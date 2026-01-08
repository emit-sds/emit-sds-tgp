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
import os, filecmp
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
    parser.add_argument('--plume_buffer_px', type=int, default=100, help='number of pixels to buffer plume cutouts by')
    parser.add_argument('--write_dcid_tifs', action='store_true', help='write out dcid level tifs for debugging')
    parser.add_argument('--n_cores', type=int, default=1, help='number of CPUs to use')
    parser.add_argument('--num_dcids', type=int, default=-1, help='number of DCIDs to process, -1 for all')
    parser.add_argument('--specific_pid', type=str, default=None, help='Run this and only this plume ID (for debugging)')
    parser.add_argument('--sync_results', action='store_true', help='sync results to remove server')
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
    fn = Filenames(args, create=True)
    fn.coverage_size = os.path.getsize(args.track_coverage_file)
    fn.annotation_size = None
    coverage = None


    # Loop only serves to rerun same instances over and over
    for run in range(max_runs):
        logging.debug('Loading Data')
        ######## Step 1 ###########

        # Update annotations file load
        utils.print_and_call(f'curl "https://popo.jpl.nasa.gov/mmgis/API/files/getfile" -H "Authorization:Bearer {args.key}" --data-raw "id={args.id}" > {fn.annotation_file_raw} 2>/dev/null')
        
        if fn.annotation_size is not None and os.path.getsize(fn.annotation_file_raw) == fn.annotation_size:
            time.sleep(10)
            continue
        fn.annotation_size = os.path.getsize(fn.annotation_file_raw)

        manual_annotations = json.load(open(fn.annotation_file_raw,'r'))['body']['geojson']
        for _feat in range(len(manual_annotations['features'])):
            manual_annotations['features'][_feat]['properties']['Plume ID'] = manual_annotations['features'][_feat]['properties'].pop('name')

        manual_annotations_previous = None
        if os.path.isfile(fn.previous_annotation_file):
            manual_annotations_previous = json.load(open(fn.previous_annotation_file,'r'))

        logging.debug('Load pre-existing public-facing outputs (combined)')
        outdict = {"crs": {"properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84" }, "type": "name"},"features":[],"name":"methane_metadata","type":"FeatureCollection" }
        if os.path.isfile(fn.output_json_internal):
            outdict = json.load(open(fn.output_json_internal,'r'))

        # If a plume didn't make it into the _internal.json, then it wasn't processed yet - remove it from the previous annotations to ensure rerun
        if len(outdict['features']) > 0 and manual_annotations_previous is not None:
            processed_plume_ids = pd.json_normalize(outdict['features'])['properties.Plume ID'].tolist()
            manual_annotations_previous['features'] = [feat for feat in manual_annotations_previous['features'] if feat['properties']['Plume ID'] in processed_plume_ids]

        # Reload coverage file if needed
        if coverage is None or os.path.getsize(args.track_coverage_file) != fn.coverage_size:
            fn.coverage_size = os.path.getsize(args.track_coverage_file)
            coverage = json.load(open(args.track_coverage_file,'r'))

      
        logging.info(f'Total plumes: {len(manual_annotations["features"])}')
        logging.debug('Run Spatial-Temporal Intersection to find FIDs')
        if args.specific_pid is None:
            manual_annotations, new_plumes = filter.add_fids(manual_annotations, coverage, manual_annotations_previous)
        else:
            # This is a special case primarily for development
            logging.info(f'Filtering to only {args.specific_pid}')
            manual_annotations['features'] = [feat for feat in manual_annotations['features'] if feat['properties']['Plume ID'] == args.specific_pid]
            manual_annotations, new_plumes = filter.add_fids(manual_annotations, coverage, None)

        logging.info(f'New plumes, post fid filter: {len(new_plumes)}')
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
        del coverage, coverage_df


        ######## Step 2 ###########
        unique_fids = np.unique([sublist for x in new_plumes for sublist in manual_annotations['features'][x]['properties']['fids']])
        unique_orbits = np.unique([manual_annotations['features'][x]['properties']['orbit'] for x in new_plumes]).tolist()
        unique_dcids = np.unique([manual_annotations['features'][x]['properties']['dcid'] for x in new_plumes]).tolist()

        if args.num_dcids > 0:
            unique_dcids = unique_dcids[0:args.num_dcids]

        # Dump out the udpated manual annotations set, so it holds FIDs / orbits for next round
        logging.debug('Dump out the updated manual annotations set, overwrting original')
        plume_io.write_geojson_linebyline(fn.annotation_file, manual_annotations)

        for feat in manual_annotations['features']:
            if 'dcid' not in feat['properties'].keys():
                logging.error('Feature missing DCID: {feat}')

        updated_windspeeds = []
        if args.n_cores > 1:
            import ray

            @ray.remote(num_cpus=1)
            def process_dcid_ray(dcid, manual_annotations, new_plumes, fn, args):
                return process_dcid(dcid, manual_annotations, new_plumes, fn, args)
            
            ray.init(num_cpus=args.n_cores, object_store_memory=10*1024**3, _temp_dir='/local/ray', include_dashboard=False)

            # Steps 3 and 4 in parallel
            jobs = [process_dcid_ray.remote(dcid, manual_annotations, new_plumes, fn, args) for dcid in unique_dcids]
            results =  ray.get(jobs)

            # Step 5
            for res in results:
                updated_plumes_point, updated_plumes_poly, ws_update = res
                outdict = update_features(outdict, updated_plumes_poly, is_point=False, add_imgs=True, fn=fn, gtype=args.type)
                outdict = update_features(outdict, updated_plumes_point, is_point=True, add_imgs=True, fn=fn, gtype=args.type)
                updated_windspeeds.extend(ws_update)

            # Final write
            plume_io.write_geojson_linebyline(fn.output_json_internal, outdict)
            plume_io.write_external_geojson(fn.output_json_internal, fn.output_json_external)
            utils.print_and_call(f'cp {fn.annotation_file} {fn.previous_annotation_file}')
            ray.shutdown()

        else:
            # Preserve the ray-free mode, mostly for debugging

            # Now step through each DCID, get the relevant FIDs, mosaic, and cut out each plume
            for dcid in unique_dcids:

                # Step 3 and 4
                updated_plumes_point, updated_plumes_poly, ws_update = process_dcid(dcid, manual_annotations, new_plumes, fn, args)
                updated_windspeeds.extend(ws_update)

                # Step 5
                outdict = update_features(outdict, updated_plumes_poly, is_point=False, add_imgs=True, fn=fn, gtype=args.type)
                outdict = update_features(outdict, updated_plumes_point, is_point=True, add_imgs=True, fn=fn, gtype=args.type)

                utils.print_and_call(f'cp {fn.annotation_file} {fn.previous_annotation_file}')
                plume_io.write_geojson_linebyline(fn.output_json_internal, outdict) # Final write
                plume_io.write_external_geojson(fn.output_json_internal, fn.output_json_external)
        
        updated_windspeeds = [x for x in updated_windspeeds if x is not None]
        if os.path.isfile(fn.working_windspeed_csv):
            wd_df = pd.read_csv(fn.working_windspeed_csv)
            ws_df = pd.concat([wd_df, pd.DataFrame(updated_windspeeds)], ignore_index=True)
        else:
            ws_df = pd.DataFrame(updated_windspeeds)
        ws_df.drop_duplicates(subset=['plume_id'], keep='last', inplace=True)
        ws_df.to_csv(fn.working_windspeed_csv, index=False)

        fn.daac_sync()

        # Sync
        if args.sync_results:
            utils.print_and_call(f'rclone copy  {fn.output_json_internal} redhat:/data/emit/mmgis/coverage/')
            utils.print_and_call(f'rclone copy  {fn.output_json_external} redhat:/data/emit/mmgis/coverage/')

    fn.daac_sync()
        


class Filenames:
    def __init__(self, args, create=False):

        self.annotation_file_raw = os.path.join(args.out_dir, "manual_annotation_raw.json") # Straight from MMGIS
        self.annotation_file = os.path.join(args.out_dir, "manual_annotation.json") # MMGIS + FIDs + Orbits 
        self.previous_annotation_file = os.path.join(args.out_dir, "previous_manual_annotation.json") # Last editted version of annotation file

        # Output / intput jsons are identical, internal just includes additional information
        self.output_json_external = os.path.join(args.out_dir, 'combined_plume_metadata_external.json') # Output (public facing) metadata file
        self.output_json_internal = os.path.join(args.out_dir, 'combined_plume_metadata_internal.json') # Output (internal facing) metadata file
        self.daac_dir = os.path.join(args.out_dir, 'daac') # Ready for DAAC sync
        self.delivery_dir = os.path.join(args.out_dir, 'delivery') # Delivery file directory
        self.quant_dir = os.path.join(args.out_dir, 'quantification') # Quantification working directory
        self.proc_dir = os.path.join(args.out_dir, 'processing') # Processing working directory
        self.working_windspeed_csv = os.path.join(args.out_dir, 'working_windspeed_estimates.csv') # Quantification windspeed working file

        if create:
            os.makedirs(self.delivery_dir, exist_ok=True)
            os.makedirs(self.quant_dir, exist_ok=True)
            os.makedirs(self.proc_dir, exist_ok=True)
            os.makedirs(self.daac_dir, exist_ok=True)
        
    @staticmethod
    def plume_delivery_basename(outdir, feat):
        return os.path.join(outdir, feat['properties']['Scene FIDs'][0][4:12], feat['properties']['Scene FIDs'][0] + '_' + feat['properties']['Plume ID'])

    @staticmethod
    def plume_working_basename(outdir, feat):
        return os.path.join(outdir, feat['properties']['Plume ID'])
 
    def feature_filenames(self, feat):
        outbase = self.plume_working_basename(self.proc_dir, feat)
        outmask_finepoly_file = os.path.join(outbase + '_finepolygon.json')
        outmask_poly_file = os.path.join(outbase + '_polygon.json')
        outmask_ort_file = os.path.join(outbase + '_mask_ort.tif')
        return outmask_finepoly_file, outmask_poly_file, outmask_ort_file

    def quantification_filenames(self, poly_plume, gtype):
        base = self.plume_working_basename(self.quant_dir, poly_plume)

        raster_file = base + '_unmasked_' + gtype.upper() + '.tif'
        unc_file = base + '_unmasked_' + gtype.upper() + '_unc.tif'
        sns_file = base + '_unmasked_' + gtype.upper() + '_sns.tif'

        return raster_file, unc_file, sns_file
    
    def delivery_filenames(self, poly_plume, gtype, daac_version=False):
        if daac_version:
            delivery_base = self.plume_delivery_basename(self.daac_dir, poly_plume)
        else:
            delivery_base = self.plume_delivery_basename(self.delivery_dir, poly_plume)

        os.makedirs(os.path.dirname(delivery_base), exist_ok=True)
        delivery_raster_file = os.path.join(delivery_base + '.tif')
        delivery_ql_file = os.path.join(delivery_base + '.png')
        delivery_json_file = os.path.join(delivery_base + '.json')
        delivery_uncert_file = delivery_raster_file.replace(gtype.upper(), gtype.upper() + '_UNC')
        delivery_sens_file = delivery_raster_file.replace(gtype.upper(), gtype.upper() + '_SNS')

        return delivery_raster_file, delivery_ql_file, delivery_json_file, delivery_uncert_file, delivery_sens_file
    
    def mmgis_q_img_filename_dict(self, plume, gtype):

        basedir = f'Layers/mosaics/images/{gtype}_q/'
        img_dicts = [{
                'name': 'Quantification Mask',
                'url': basedir +  plume['properties']['Plume ID'] + '_ql.png',
                'type': 'image'
                },
                {
                'name': 'Quantification Sweep',
                'url': basedir +  plume['properties']['Plume ID'] + '_sweep.pdf',
                'type': 'document'
            }]
        return img_dicts

    def daac_sync(self):

        output_plumes = [x for x in json.load(open(self.output_json_external, 'r'))['features'] if x['geometry']['type'] == 'Polygon']
        for plume in output_plumes:
            delivery_raster, delivery_ql, delivery_json, delivery_uncert, delivery_sens = self.delivery_filenames(plume, 'point', daac_version=False)
            daac_raster, daac_ql, daac_json, daac_uncert, daac_sens = self.delivery_filenames(plume, 'point', daac_version=True)
            for de, da in zip([delivery_raster, delivery_ql, delivery_json],
                               [daac_raster, daac_ql, daac_json]):

                if os.path.isfile(da):
                    # If it's there and the same, skip safely
                    if filecmp.cmp(de, da, shallow=False):
                        continue
                    else:
                        # If it's there and different, raise the alarm
                        logging.warning(f'DAAC sync - file {da} exists and is different from delivery - we cannot delivery this file in the same version')


                else:
                    # If it's not there yet, copy
                    utils.print_and_call(f'cp {de} {da}')
                    continue







def process_dcid(dcid, manual_annotations, new_plumes, fn, args):

    logging.info(f'Processing DCID {dcid}...')
    #plumes_idx_in_dcid = [x for x in range(len(manual_annotations['features'])) if manual_annotations['features'][x]['properties']['dcid'] == dcid]
    plumes_idx_in_dcid = [x for x in new_plumes if manual_annotations['features'][x]['properties']['dcid'] == dcid]
    fids_in_dcid = np.unique([sublist for x in plumes_idx_in_dcid for sublist in manual_annotations['features'][x]['properties']['fids']])
    logging.info(f'...found {len(plumes_idx_in_dcid)} plumes to update and {len(fids_in_dcid)} FIDs')

    ort_dat_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type) for fid in fids_in_dcid]
    ort_sens_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type, data_value='sens') for fid in fids_in_dcid]
    unc_dat_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type, data_value='uncert') for fid in fids_in_dcid]

    dcid_ort_vrt_file = os.path.join(fn.proc_dir, f'dcid_{dcid}_mf_ort.vrt')
    dcid_ort_unc_vrt_file = os.path.join(fn.proc_dir, f'dcid_{dcid}_unc_ort.vrt')
    dcid_ort_sns_vrt_file = os.path.join(fn.proc_dir, f'dcid_{dcid}_sns_ort.vrt')

    dcid_ort_tif_file = os.path.join(fn.proc_dir, f'dcid_{dcid}_mf_ort.tif')
    dcid_ort_unc_tif_file = os.path.join(fn.proc_dir, f'dcid_{dcid}_unc_ort.tif')
    dcid_ort_sns_tif_file = os.path.join(fn.proc_dir, f'dcid_{dcid}_sns_ort.tif')
    utils.print_and_call(f'gdalbuildvrt {dcid_ort_vrt_file} {" ".join(ort_dat_files)} --quiet')
    if args.write_dcid_tifs:
        utils.print_and_call(f'gdal_translate {dcid_ort_vrt_file} {dcid_ort_tif_file} -co COMPRESS=LZW --quiet')
    else:
        dcid_ort_tif_file = dcid_ort_vrt_file  # just use the VRT directly to save space/time

    ort_ds = gdal.Open(dcid_ort_tif_file)
    trans = ort_ds.GetGeoTransform()

    # Only create sns and unc if needed

    utils.print_and_call(f'gdalbuildvrt {dcid_ort_unc_vrt_file} {" ".join(unc_dat_files)} --quiet')
    utils.print_and_call(f'gdalbuildvrt {dcid_ort_sns_vrt_file} {" ".join(ort_sens_files)} --quiet')
    if args.write_dcid_tifs:
        utils.print_and_call(f'gdal_translate {dcid_ort_unc_vrt_file} {dcid_ort_unc_tif_file} -co COMPRESS=LZW --quiet')
        utils.print_and_call(f'gdal_translate {dcid_ort_sns_vrt_file} {dcid_ort_sns_tif_file} -co COMPRESS=LZW --quiet')
    else:
        dcid_ort_unc_tif_file = dcid_ort_unc_vrt_file
        dcid_ort_sns_tif_file = dcid_ort_sns_vrt_file
    unc_ds = gdal.Open(dcid_ort_unc_tif_file)
    sns_ds = gdal.Open(dcid_ort_sns_tif_file)

    # Calculate pixel size for DCID
    proj_ds = gdal.Warp('', dcid_ort_tif_file, dstSRS='EPSG:3857', format='VRT')
    transform_3857 = proj_ds.GetGeoTransform()
    xsize_m = transform_3857[1]
    ysize_m = transform_3857[5]
    del proj_ds


    ######## Step 3 ##############
    # Use the manual plumes to come up with a new set of plume masks and labels
    updated_plumes_poly, updated_plumes_point, updated_windspeed = [], [], []
    for newp in plumes_idx_in_dcid:
        feat = manual_annotations['features'][newp]
        logging.info(f'Building output plume {feat["properties"]["Plume ID"]}')

        deliver_emissions = feat['properties']['Simple IME Valid'] == 'Yes'

        rawspace_coords = plume_io.rawspace_coordinate_conversion([], feat['geometry']['coordinates'][0], trans, ortho=True)

        datshape = (ort_ds.RasterYSize, ort_ds.RasterXSize)
        window, newp_trans, local_coords = plume_io.get_window(rawspace_coords, trans, datshape, args.plume_buffer_px)
        if window is None:
            logging.warning(f'Plume {feat["properties"]["Plume ID"]} has invalid window, skipping')
            continue

        cut_plume_data = ort_ds.ReadAsArray(window[0], window[1], window[2], window[3]).squeeze()
        cut_uncdat = unc_ds.ReadAsArray(window[0], window[1], window[2], window[3]).squeeze()
        cut_snsdat = sns_ds.ReadAsArray(window[0], window[1], window[2], window[3]).squeeze()

        manual_mask = rasterize(shapes=[Polygon(local_coords)], out_shape=(cut_plume_data.shape[0],cut_plume_data.shape[1]), dtype=np.uint8) # numpy binary mask for manual IDs

        plumestyle = 'classic'
        if 'Delineation Mode' in feat['properties'].keys():
            plumestyle = feat['properties']['Delineation Mode']

        loc_fid_mask = None # change name to dcid_plume_mask
        if plumestyle == 'classic':
            loc_fid_mask = utils.plume_mask_threshold(cut_plume_data.copy(), manual_mask, style=args.type)
        elif plumestyle == 'manual':    
            loc_fid_mask = manual_mask.astype(bool)
        

        ############  Step 4 ###########
        outmask_finepoly_file, outmask_poly_file, outmask_ort_file = fn.feature_filenames(feat) 

        # Write mask file and save tif reference
        #write_output_file(newp_trans, ort_ds.GetProjection(), cut_plume_mask, outmask_ort_file)
        #write_output_file(newp_trans, ort_ds.GetProjection(), cut_plume_data, outmask_ort_file)
        plume_io.write_cog(outmask_ort_file, loc_fid_mask.reshape((loc_fid_mask.shape[0], loc_fid_mask.shape[1],1)).astype(np.uint8), 
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



        delivery_raster_file, delivery_ql_file, delivery_json_file, delivery_uncert_file, delivery_sens_file = fn.delivery_filenames(poly_plume, args.type)

        # Write delivery files
        meta = plume_io.get_metadata(poly_plume, plume_io.global_metadata(data_version=args.data_version))
        plume_io.write_cog(delivery_raster_file, cut_plume_data.astype(np.float32), newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta, mask=loc_fid_mask)
        plume_io.write_cog(delivery_uncert_file, cut_uncdat.astype(np.float32), newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta, mask=loc_fid_mask)
        plume_io.write_cog(delivery_sens_file, cut_snsdat.astype(np.float32), newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta, mask=loc_fid_mask)
        plume_io.write_color_quicklook(cut_plume_data, delivery_ql_file, inmask=loc_fid_mask)

        # Write unmasked version of delivery files for quantification (mainly for plotting)
        quant_raster_file, quant_uncert_file, quant_sens_file = fn.quantification_filenames(poly_plume, args.type)
        plume_io.write_cog(quant_raster_file, cut_plume_data.astype(np.float32), newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)
        plume_io.write_cog(quant_uncert_file, cut_uncdat.astype(np.float32), newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)
        plume_io.write_cog(quant_sens_file, cut_snsdat.astype(np.float32), newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)

        # Compute Emissions
        emissions_info, windspeed_info = compute_Q_and_uncertainty_utils.single_plume_emissions(
            feat,
            poly_plume,
            fn.quant_dir,
            fn.proc_dir,
            quant_raster_file,
            quant_sens_file,
            quant_uncert_file,
            fn.annotation_file,
            working_windspeed_csv=fn.working_windspeed_csv,
            overrule_simple_ime_flag=True, # we want to run the calc no matter what - we'll discard later per metadata
        )


        poly_plume['properties'].update(emissions_info)
        point_plume['properties'].update(emissions_info)

        # For the delivery file, if not flagged for emissions delivery, don't include
        plume_io.write_delivery_json(delivery_json_file, poly_plume, meta['DAAC Scene Names'], deliver_emissions)

        # Now save for output / archive jsons
        updated_plumes_point.append(point_plume)
        updated_plumes_poly.append(poly_plume)
        if windspeed_info is not None:
            updated_windspeed.append(windspeed_info)


    return updated_plumes_point, updated_plumes_poly, updated_windspeed


def update_features(existing: dict, new_features: List, is_point: bool, add_imgs: dict=None, fn: Filenames=None, gtype: str = 'ch4') -> dict:

    if int(add_imgs is not None) + int(fn is not None) == 1:
        logging.warning('Both add_imgs and fn must be provided, or neither; skipping image addition')

    ########## Step 5 ##########
    for plm in new_features:
        if is_point:
            existing_match_index = [_x for _x, x in enumerate(existing['features']) if plm['properties']['Plume ID'] == x['properties']['Plume ID'] and x['geometry']['type'] == 'Point']
        else:
            existing_match_index = [_x for _x, x in enumerate(existing['features']) if plm['properties']['Plume ID'] == x['properties']['Plume ID'] and x['geometry']['type'] != 'Point']
            
        if len(existing_match_index) > 2:
            logging.warning("HELP! Too many matching indices")

        if add_imgs is not None and fn is not None:
            plm['properties']['images'] = fn.mmgis_q_img_filename_dict(plm, gtype.lower())

        if len(existing_match_index) > 0:
            existing['features'][existing_match_index[0]] = plm
        else:
            existing['features'].append(plm)
        
        
    return existing
 





 
if __name__ == '__main__':
    main()



