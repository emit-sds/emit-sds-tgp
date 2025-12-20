#! /usr/bin/env python
#
#  Copyright 2023 California Institute of Technology
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
# Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov

'''
Major Steps, at each iteration of the loop:

1. Download the latest plume annotations, and intersect them with coverage to add FIDS, Orbits, etc.  ID new plumes during this
2. Find all DCIDs that have new plumes, and mosaic each 
3. Step through each plume in the DCID, cut it out, calculate stats and metadata
4. Prep each plume for delivery, writing output files
5. Merge all plumes into single output metadata, merge all plumes within a DCID to a single COG / tiled dataset

Note that previously, the necessity of the background calculation required an extra iteration through each DCID in 3, which was costly.
Removing the need for the background calcluation eliminates this step, and it has been cut out of the code for speed.
'''


import argparse
import os
import numpy as np
import json
import time
from shapely.geometry import Polygon
from osgeo import gdal
from rasterio.features import rasterize
import logging
from emit_main.workflow.workflow_manager import WorkflowManager
import pandas as pd
from typing import List

from annotate import plume_io, filter, utils
from quantification import compute_flux, windspeed, compute_Q_and_uncertainty_utils

gdal.osr.UseExceptions()


def get_sds_cog(fid, enh_version, dtype='ch4', data_value=''):
    date = fid[4:12]
    path=f'/store/emit/ops/data/acquisitions/{date}/{fid.split("_")[0]}/ghg/{dtype}/{fid.split("_")[0]}*_ghg_ort{data_value}{dtype}_b0106_{enh_version}.tif'
    return path

def plume_delivery_basename(outdir, feat):
    return os.path.join(outdir, feat['properties']['Scene FIDs'][0][4:12], feat['properties']['Scene FIDs'][0] + '_' + feat['properties']['Plume ID'])

def plume_working_basename(outdir, feat):
    return os.path.join(outdir, feat['properties']['Plume ID'])
    


def main(input_args=None):
    parser = argparse.ArgumentParser(description="merge jsons")
    parser.add_argument('key', type=str,  metavar='INPUT_DIR', help='input directory')   
    parser.add_argument('id', type=str,  metavar='INPUT_DIR', help='input directory')   
    parser.add_argument('out_dir', type=str,  metavar='OUTPUT_DIR', help='output directory')   
    parser.add_argument('data_version', type=str)
    parser.add_argument('--enh_data_version', type=str, default='v02')
    parser.add_argument('--type', type=str,  choices=['ch4','co2'], default='ch4')   
    parser.add_argument('--database_config', type=str,  default='/store/emit/ops/repos/emit-main/emit_main/config/ops_sds_config.json')   
    parser.add_argument('--loglevel', type=str, default='DEBUG', help='logging verbosity')    
    parser.add_argument('--logfile', type=str, default=None, help='output file to write log to')    
    parser.add_argument('--continuous', action='store_true', help='run continuously')    
    parser.add_argument('--track_coverage_file', default='/store/brodrick/emit/emit-visuals/track_coverage_pub.json')
    parser.add_argument('--plume_buffer_px', type=int, default=30, help='number of pixels to buffer plume cutouts by')
    parser.add_argument('--write_dcid_tifs', action='store_true', help='write out dcid level tifs for debugging')
    args = parser.parse_args(input_args)

    logging.basicConfig(format='%(levelname)s:%(asctime)s ||| %(message)s', level=args.loglevel,
                        filename=args.logfile, datefmt='%Y-%m-%d,%H:%M:%S')

    np.random.seed(13)
    max_runs = 1
    if args.continuous:
        max_runs = int(1e15)

    try:
        database = WorkflowManager(config_path=args.database_config).database_manager
    except:
        raise AttributeError('Could not open databse - check args.database_config')

    # Global File names
    annotation_file_raw = os.path.join(args.out_dir, "manual_annotation_raw.json") # Straight from MMGIS
    annotation_file = os.path.join(args.out_dir, "manual_annotation.json") # MMGIS + FIDs + Orbits 
    previous_annotation_file = os.path.join(args.out_dir, "previous_manual_annotation.json") # Last editted version of annotation file
    output_json = os.path.join(args.out_dir, 'combined_plume_metadata.json') # Output (public facing) metadata file
    delivery_dir = os.path.join(args.out_dir, 'delivery') # Delivery file directory
    quant_dir = os.path.join(args.out_dir, 'quantification') # Quantification working directory
    proc_dir = os.path.join(args.out_dir, 'processing') # Processing working directory
    working_windspeed_csv = os.path.join(args.out_dir, 'working_windspeed_estimates_0000.csv') # Quantification windspeed working file
    os.makedirs(delivery_dir, exist_ok=True)
    os.makedirs(quant_dir, exist_ok=True)
    os.makedirs(proc_dir, exist_ok=True)


    # Loop only serves to rerun same instances over and over
    for run in range(max_runs):
        logging.debug('Loading Data')
        ######## Step 1 ###########

        # Update annotations file load
        utils.print_and_call(f'curl "https://popo.jpl.nasa.gov/mmgis/API/files/getfile" -H "Authorization:Bearer {args.key}" --data-raw "id={args.id}" > {annotation_file_raw}')
        manual_annotations = json.load(open(annotation_file_raw,'r'))['body']['geojson']
        for _feat in range(len(manual_annotations['features'])):
            manual_annotations['features'][_feat]['properties']['Plume ID'] = manual_annotations['features'][_feat]['properties'].pop('name')
        manual_annotations_previous = None
        if os.path.isfile(previous_annotation_file):
            manual_annotations_previous = json.load(open(previous_annotation_file,'r'))


        # Reload coverage file load
        coverage = json.load(open(args.track_coverage_file,'r'))

        logging.debug('Load pre-existing public-facing outputs (combined)')
        outdict = {"crs": {"properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84" }, "type": "name"},"features":[],"name":"methane_metadata","type":"FeatureCollection" }
        if os.path.isfile(output_json):
            outdict = json.load(open(output_json,'r'))

        logging.debug('Run Spatial-Temporal Intersection to find FIDs')
        manual_annotations, new_plumes = filter.add_fids(manual_annotations, coverage, manual_annotations_previous)
        logging.info(f'New plume, post fid: {len(new_plumes)}')
        manual_annotations_df = pd.json_normalize(manual_annotations['features'])
        most_recent_plume_create = pd.to_datetime(manual_annotations_df['properties.Time Created']).max()
        last_time_range = pd.to_datetime(manual_annotations_df['properties.Time Range End']).dt.tz_localize('UTC').max()


        logging.debug('Pull Summary Stats')
        coverage_df = pd.json_normalize(coverage['features'])
        r0_review_needed = np.sum(pd.to_datetime(coverage_df['properties.end_time']) > last_time_range)

        double_approved_count = np.sum([ np.all([x['properties'][k] for k in ['R1 - Reviewed', 'R1 - VISIONS', 'R2 - Reviewed', 'R2 - VISIONS']]) for x in manual_annotations['features']])
        populated_emission_count = np.sum([ np.all([x['properties'][k] for k in ['R1 - Reviewed', 'R1 - VISIONS', 'R2 - Reviewed', 'R2 - VISIONS']]) and x['properties']['Simple IME Valid'] != 'NA' for x in manual_annotations['features']])
        reporting_emission_count = np.sum([ np.all([x['properties'][k] for k in ['R1 - Reviewed', 'R1 - VISIONS', 'R2 - Reviewed', 'R2 - VISIONS']]) and x['properties']['Simple IME Valid'] == 'Yes' for x in manual_annotations['features']])
        r1_review_count = np.sum([ x['properties']['R1 - Reviewed'] is False for x in manual_annotations['features']])
        r2_review_count = np.sum([ np.all([x['properties'][k] for k in ['R1 - Reviewed', 'R1 - VISIONS']]) and not x['properties']['R2 - Reviewed'] for x in manual_annotations['features']])
        logging.info(f'Plume Complexes Approved for VISIONS: {double_approved_count}')
        logging.info(f'Approved Plume Complexes with Simple IME Evaluated: {populated_emission_count}')
        logging.info(f'Plume Complexes with Emission Estimate: {reporting_emission_count}')
        logging.info(f'R1 Review Deck: {r1_review_count}')
        logging.info(f'R2 Review Deck: {r2_review_count}')
        logging.info(f'Most Recent R0 Plume Created: {most_recent_plume_create}')
        logging.info(f'Most Recent R0 Plume Scene: {last_time_range}')
        logging.info(f'Most Recent Scene: {coverage["features"][-1]["properties"]["end_time"]}')
        logging.info(f'Inferred scenes needing R0 Review: {r0_review_needed}')
        # If there's nothing new, sleep and retry
        if len(new_plumes) == 0:
            time.sleep(10)
            continue

        logging.debug('Querry the EMIT DB to add Orbits and DCIDs')
        manual_annotations, new_plumes     = filter.add_orbits(manual_annotations, new_plumes, database)
        logging.info(f'New plume, post orbit: {len(new_plumes)}')



        ######## Step 2 ###########
        unique_fids = np.unique([sublist for x in new_plumes for sublist in manual_annotations['features'][x]['properties']['fids']])
        unique_orbits = np.unique([manual_annotations['features'][x]['properties']['orbit'] for x in new_plumes]).tolist()
        unique_dcids = np.unique([manual_annotations['features'][x]['properties']['dcid'] for x in new_plumes]).tolist()

        # Dump out the udpated manual annotations set, so it holds FIDs / orbits for next round
        logging.debug('Dump out the updated manual annotations set, overwrting original')
        plume_io.write_geojson_linebyline(annotation_file, manual_annotations)
        #with open(annotation_file, 'w') as fout:
        #    fout.write(json.dumps(manual_annotations, cls=plume_io.SerialEncoder)) 

        for feat in manual_annotations['features']:
            if 'dcid' not in feat['properties'].keys():
                logging.error('Feature missing DCID: {feat}')


        # Now step through each DCID, get the relevant FIDs, mosaic, and cut out each plume
        for dcid in unique_dcids:

            logging.info(f'Processing DCID {dcid}...')
            #plumes_idx_in_dcid = [x for x in range(len(manual_annotations['features'])) if manual_annotations['features'][x]['properties']['dcid'] == dcid]
            plumes_idx_in_dcid = [x for x in new_plumes if manual_annotations['features'][x]['properties']['dcid'] == dcid]
            fids_in_dcid = np.unique([sublist for x in plumes_idx_in_dcid for sublist in manual_annotations['features'][x]['properties']['fids']])
            logging.info(f'...found {len(plumes_idx_in_dcid)} plumes to update and {len(fids_in_dcid)} FIDs')

            ort_dat_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type) for fid in fids_in_dcid]
            ort_sens_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type, data_value='sens') for fid in fids_in_dcid]
            unc_dat_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type, data_value='uncert') for fid in fids_in_dcid]

            dcid_ort_vrt_file = os.path.join(proc_dir, f'dcid_{dcid}_mf_ort.vrt')
            dcid_ort_unc_vrt_file = os.path.join(proc_dir, f'dcid_{dcid}_unc_ort.vrt')
            dcid_ort_sns_vrt_file = os.path.join(proc_dir, f'dcid_{dcid}_sns_ort.vrt')

            dcid_ort_tif_file = os.path.join(proc_dir, f'dcid_{dcid}_mf_ort.tif')
            dcid_ort_unc_tif_file = os.path.join(proc_dir, f'dcid_{dcid}_unc_ort.tif')
            dcid_ort_sns_tif_file = os.path.join(proc_dir, f'dcid_{dcid}_sns_ort.tif')
            utils.print_and_call(f'gdalbuildvrt {dcid_ort_vrt_file} {" ".join(ort_dat_files)} --quiet')
            if args.write_dcid_tifs:
                utils.print_and_call(f'gdal_translate {dcid_ort_vrt_file} {dcid_ort_tif_file} -co COMPRESS=LZW --quiet')
            else:
                dcid_ort_tif_file = dcid_ort_vrt_file  # just use the VRT directly to save space/time

            ort_ds = gdal.Open(dcid_ort_tif_file)
            ortdat = ort_ds.ReadAsArray().squeeze()
            trans = ort_ds.GetGeoTransform()

            # Only create sns and unc if needed
            estimate_simple_ime = np.sum([manual_annotations['features'][x]['properties']['Simple IME Valid'] == 'Yes' for x in plumes_idx_in_dcid]) > 0
            snsdat, uncdat = None, None
            if estimate_simple_ime:
                utils.print_and_call(f'gdalbuildvrt {dcid_ort_unc_vrt_file} {" ".join(unc_dat_files)} --quiet')
                utils.print_and_call(f'gdalbuildvrt {dcid_ort_sns_vrt_file} {" ".join(ort_sens_files)} --quiet')
                if args.write_dcid_tifs:
                    utils.print_and_call(f'gdal_translate {dcid_ort_unc_vrt_file} {dcid_ort_unc_tif_file} -co COMPRESS=LZW --quiet')
                    utils.print_and_call(f'gdal_translate {dcid_ort_sns_vrt_file} {dcid_ort_sns_tif_file} -co COMPRESS=LZW --quiet')
                else:
                    dcid_ort_unc_tif_file = dcid_ort_unc_vrt_file
                    dcid_ort_sns_tif_file = dcid_ort_sns_vrt_file
                uncdat = gdal.Open(dcid_ort_unc_tif_file).ReadAsArray().squeeze()
                snsdat = gdal.Open(dcid_ort_sns_tif_file).ReadAsArray().squeeze()

            # Calculate pixel size for DCID
            proj_ds = gdal.Warp('', dcid_ort_tif_file, dstSRS='EPSG:3857', format='VRT')
            transform_3857 = proj_ds.GetGeoTransform()
            xsize_m = transform_3857[1]
            ysize_m = transform_3857[5]
            del proj_ds




            ######## Step 3 ##############
            # Use the manual plumes to come up with a new set of plume masks and labels
            combined_mask = None
            updated_plumes_poly, updated_plumes_point = [], [] 
            for newp in plumes_idx_in_dcid:
                feat = manual_annotations['features'][newp]
                logging.info(f'Building output plume {feat["properties"]["Plume ID"]}')

                rawspace_coords = plume_io.rawspace_coordinate_conversion([], feat['geometry']['coordinates'][0], trans, ortho=True)
                manual_mask = rasterize(shapes=[Polygon(rawspace_coords)], out_shape=(ortdat.shape[0],ortdat.shape[1])) # numpy binary mask for manual IDs

                plumestyle = 'classic'
                if 'Delineation Mode' in feat['properties'].keys():
                    plumestyle = feat['properties']['Delineation Mode']

                loc_fid_mask = None # change name to dcid_plume_mask
                if plumestyle == 'classic':
                    loc_fid_mask = utils.plume_mask_threshold(ortdat.copy(), manual_mask, style=args.type)
                elif plumestyle == 'manual':    
                    loc_fid_mask = manual_mask.astype(bool)
                
                if combined_mask is None:
                    combined_mask = loc_fid_mask.copy()
                else:
                    combined_mask = np.logical_or(combined_mask, loc_fid_mask)

                cut_plume_mask, newp_trans, s = plume_io.trim_plume(loc_fid_mask, trans, buffer=args.plume_buffer_px)
                cut_plume_data, _, s = plume_io.trim_plume(ortdat, trans, badmask=loc_fid_mask == 0, nodata_value=-9999, buffer=args.plume_buffer_px)
                if s is False:
                    logging.warning(f'Plume {feat["properties"]["Plume ID"]} has no valid pixels after trimming, skipping plume')
                    continue

                if estimate_simple_ime:
                    cut_uncdat, _, s = plume_io.trim_plume(uncdat, trans, badmask=loc_fid_mask == 0, nodata_value=-9999, buffer=args.plume_buffer_px)
                    cut_snsdat, _, s = plume_io.trim_plume(snsdat, trans, badmask=loc_fid_mask == 0, nodata_value=-9999, buffer=args.plume_buffer_px)



                ############  Step 4 ###########
                outbase = plume_working_basename(proc_dir, feat)
                outmask_finepoly_file = os.path.join(outbase + '_finepolygon.json')
                outmask_poly_file = os.path.join(outbase + '_polygon.json')
                outmask_ort_file = os.path.join(outbase + '_mask_ort.tif')

                # Write mask file and save tif reference
                #write_output_file(newp_trans, ort_ds.GetProjection(), cut_plume_mask, outmask_ort_file)
                #write_output_file(newp_trans, ort_ds.GetProjection(), cut_plume_data, outmask_ort_file)
                plume_io.write_cog(outmask_ort_file, cut_plume_mask.reshape((cut_plume_mask.shape[0], cut_plume_mask.shape[1],1)).astype(np.uint8), 
                          newp_trans, ort_ds.GetProjection(), nodata_value=0)
                
                if os.path.isfile(outmask_poly_file):
                    os.remove(outmask_poly_file)
                if os.path.isfile(outmask_finepoly_file):
                    os.remove(outmask_finepoly_file)
                utils.print_and_call(f'gdal_polygonize {outmask_ort_file} {outmask_finepoly_file} -f GeoJSON -mask {outmask_ort_file} -8 -quiet')
                utils.print_and_call(f'ogr2ogr {outmask_poly_file} {outmask_finepoly_file} -f GeoJSON -lco RFC7946=YES -simplify {trans[1]} --quiet')

                # Read back in the polygon we just wrote
                plume_to_add = json.load(open(outmask_poly_file))['features']
                if len(plume_to_add) > 1:
                    logging.warning(f'ACK - multiple polygons from one Plume ID in file {outmask_poly_file}')
                plume_to_add[0]['geometry']['coordinates'] = [[[np.round(x[0],5), np.round(x[1],5)] for x in plume_to_add[0]['geometry']['coordinates'][0]]]

                poly_plume, point_plume = utils.build_plume_properties(feat['properties'], plume_to_add[0]['geometry'], cut_plume_data, 
                                                                       newp_trans, args.data_version, xsize_m, ysize_m)



                # These are the real delivery files
                delivery_base = plume_delivery_basename(delivery_dir, poly_plume)
                os.makedirs(os.path.dirname(delivery_base), exist_ok=True)

                delivery_raster_file = os.path.join(delivery_base + '.tif')
                delivery_ql_file = os.path.join(delivery_base + '.png')
                delivery_json_file = os.path.join(delivery_base + '.json')
                delivery_uncert_file = delivery_raster_file.replace(args.type.upper(), args.type.upper() + '_UNC')
                delivery_sens_file = delivery_raster_file.replace(args.type.upper(), args.type.upper() + '_SNS')

                meta = plume_io.get_metadata(poly_plume, plume_io.global_metadata(data_version=args.data_version))
                plume_io.write_cog(delivery_raster_file, cut_plume_data.reshape((cut_plume_data.shape[0], cut_plume_data.shape[1],1)).astype(np.float32), 
                          newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)
                plume_io.write_color_quicklook(cut_plume_data, delivery_ql_file)

                emissions_info = {
                    'Wind Speed (m/s)': 'NA',
                    'Wind Speed Std (m/s)': 'NA',
                    'Wind Speed Source': 'NA',
                    'Emissions Rate Estimate (kg/hr)': 'NA',
                    'Emissions Rate Estimate Uncertainty (kg/hr)': 'NA',
                    'Fetch Length (m)': 'NA',
                }
                if estimate_simple_ime and feat['properties']['Simple IME Valid'] == 'Yes':
                    plume_io.write_cog(delivery_uncert_file, cut_uncdat.reshape((cut_uncdat.shape[0], cut_uncdat.shape[1],1)).astype(np.float32), 
                              newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)
                    plume_io.write_cog(delivery_sens_file, cut_snsdat.reshape((cut_snsdat.shape[0], cut_snsdat.shape[1],1)).astype(np.float32), 
                              newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)


                    class flux_args:
                        def __init__(self):
                            self.minppmm = 500
                            self.mergedistm = 200
                            self.maxppmm = 3000
                            self.maxfetchm = 1000
                            self.minaream2 = 0

                            self.plot_diag = True
                            self.verbose = False
                            self.fid = poly_plume['properties']['Plume ID']

                            self.name_suffix = ''
                            self.plot_path = quant_dir

                            if feat['properties']['Psuedo-Origin'] in [None, '']:
                                self.lat = None
                                self.lng = None
                            else:
                                pso = json.loads(feat['properties']['Psuedo-Origin'])
                                self.lat = pso['coordinates'][1]
                                self.lng = pso['coordinates'][0]

                            self.cmfimgf = delivery_raster_file
                            self.sns_file = delivery_sens_file
                            self.unc_file = delivery_uncert_file
 
                            mb_coords = poly_plume['geometry']['coordinates']
                            if isinstance(mb_coords[0][0], List):
                                mb_coords = mb_coords[0] 
                            self.manual_boundary_coordinates_lon_lat = [item for sublist in mb_coords for item in sublist]

                            self.mask_mode = 'infer'
                            self.q_method = 'mask'

                            self.rgbimgf = None
                            self.labimgf = None

                            #self.log_level = args.loglevel
                            self.log_level = 'DEBUG'
                            self.logfile = None

                        unc_file = None

                    lfa = flux_args()
                    if lfa.lat is None or lfa.lng is None:
                        logging.debug(f'Plume {lfa.fid} missing Psuedo-Origin, skipping plume')
                        continue

                    with open(os.path.join(proc_dir, f'flux_args_{poly_plume["properties"]["Plume ID"]}.json'),'wb') as pf:
                        pf.write(json.dumps(lfa.__dict__, indent=2).encode('utf-8'))


                    # Compute Flux
                    #quant_res: [plume_complex, C_Q_MASK, C_Q_CC, lng, lat, fetchm, mergedistm, args.minppmm, args.maxppmm, args.minaream2, ps, C2_UNC_MASK]
                    flux_status, flux_res = compute_flux.compute_flux(lfa)
                    if flux_status != 'success':
                        logging.warning(f'Flux calculation failed for plume {lfa.fid} with status {flux_status}, skipping plume')
                        continue
                    logging.debug(f'Flux results: {flux_status} {flux_res}')

                    # Compute Windspeed
                    original_log_level = logging.getLogger().level
                    logging.getLogger().setLevel(logging.ERROR)
                    windspeed.update_EMIT_plume_list_windspeeds(working_windspeed_csv,
                                                                plume_file=annotation_file,
                                                                write_rate=1,
                                                                plume_list=[poly_plume['properties']['Plume ID']]
                                                                )
                    logging.getLogger().setLevel(original_log_level)
                    wsk = windspeed.windspeed_key_names('hrrr', 'era5', 'w10', 'm_per_s')

                    dfw = pd.read_csv(working_windspeed_csv)
                    dfw = dfw[dfw['plume_id'] == poly_plume['properties']['Plume ID']]

                    ws = dfw[wsk['primary']].fillna(dfw[wsk['secondary']])
                    ws_std = dfw[wsk['primary_std']].fillna(dfw[wsk['secondary_std']])
                    ws_source = np.where(dfw[wsk['primary']].notna(), 'hrrr', 'era5')

                    # Compute Emissions
                    Q, sigma_Q, sigma_C, sigma_w, sigma_f = \
                        compute_Q_and_uncertainty_utils.compute_Q_and_uncertainty(
                            flux_res[1], # C_Q_mask_kg_hr_mpers
                            flux_res[11], #sum_C_sqrd
                            ws, 
                            ws_std, 
                            flux_res[5], #fetch
                            0.0) #fetch_unc_frac
                    emissions_info['Wind Speed (m/s)'] = float(np.round(ws.values[0],4))
                    emissions_info['Wind Speed Std (m/s)'] = float(np.round(ws_std.values[0],4))
                    emissions_info['Wind Speed Source'] = ws_source[0].upper()
                    emissions_info['Emissions Rate Estimate (kg/hr)'] = float(np.round(Q.values[0],4))
                    emissions_info['Emissions Rate Estimate Uncertainty (kg/hr)'] = float(np.round(sigma_Q.values[0],4))

                    emissions_info['Fetch Length (m)'] = float(np.round(flux_res[5],4))
                    logging.info(f'Populated emissions info: {emissions_info}')

                poly_plume['properties'].update(emissions_info)
                point_plume['properties'].update(emissions_info)

                updated_plumes_point.append(point_plume)
                updated_plumes_poly.append(poly_plume)
                plume_io.write_json(delivery_json_file, poly_plume, meta['Source_Scenes'])

        ########## Step 5 ##########
        for plm in updated_plumes_poly:
            existing_match_index = [_x for _x, x in enumerate(outdict['features']) if plm['properties']['Plume ID'] == x['properties']['Plume ID'] and x['geometry']['type'] != 'Point']
            if len(existing_match_index) > 2:
                logging.warning("HELP! Too many matching indices")
            if len(existing_match_index) > 0:
                outdict['features'][existing_match_index[0]] = plm
            else:
                outdict['features'].append(plm)
 
        for plm in updated_plumes_point:
            existing_match_index = [_x for _x, x in enumerate(outdict['features']) if plm['properties']['Plume ID'] == x['properties']['Plume ID'] and x['geometry']['type'] == 'Point']
            if len(existing_match_index) > 2:
                logging.warning("HELP! Too many matching indices")
            if len(existing_match_index) > 0:
                outdict['features'][existing_match_index[0]] = plm
            else:
                outdict['features'].append(plm)

        plume_io.write_geojson_linebyline(output_json, outdict) # Final write


        # Sync
        utils.print_and_call(f'cp {annotation_file} {previous_annotation_file}')
        #print_and_call(f'rsync -a --info=progress2 {tile_dir}/{od_date}/ brodrick@$EMIT_SCIENCE_IP:/data/emit/mmgis/mosaics/{args.type}_plume_tiles_working/{od_date}/ --delete')
        #print_and_call(f'rsync {output_json} brodrick@$EMIT_SCIENCE_IP:/data/emit/mmgis/coverage/converted_manual_{args.type}_plumes.json')





        # Estimating the background would require going through *ALL* plumes in the DCID, not just the new ones - no longer needed, just leaving as a reminder
        #background = ortdat[np.logical_and(combined_mask == 0, ortdat != -9999)]















 
if __name__ == '__main__':
    main()



