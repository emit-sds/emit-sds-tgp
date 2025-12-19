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

import argparse
import os
import numpy as np
import json
import glob
import datetime
import time
import subprocess
from spectral.io import envi
from shapely.geometry import Polygon
import shapely
from copy import deepcopy
from osgeo import gdal
from rasterio.features import rasterize
import logging
import matplotlib.pyplot as plt
from emit_main.workflow.workflow_manager import WorkflowManager
import pandas as pd
import geopandas as gpd
from scipy.ndimage import gaussian_filter



def gauss_blur(img, sigma, preserve_nans=False):
    
    V=img.copy()
    V[np.isnan(img)]=0
    VV=gaussian_filter(V,sigma=sigma)

    W=0*img.copy()+1
    W[np.isnan(img)]=0
    WW=gaussian_filter(W,sigma=sigma)

    img_smooth=VV/WW
    if preserve_nans:
        img_smooth[np.isnan(img)] = np.nan
    
    return img_smooth

def print_and_call(cmd_str):
    logging.debug(f'Calling: {cmd_str}')
    subprocess.call(cmd_str, shell=True)

#def write_output_file(source_trans, source_proj, output_img, output_file):
#    driver = gdal.GetDriverByName('GTiff')
#    driver.Register()
#    if len(output_img.shape) == 2:
#        outDataset = driver.Create(output_file,output_img.shape[1],output_img.shape[0],1,gdal.GDT_Byte,options = ['COMPRESS=LZW'])
#        outDataset.GetRasterBand(1).WriteArray(output_img)
#    else:
#        outDataset = driver.Create(output_file,output_img.shape[1],output_img.shape[0],3,gdal.GDT_Byte,options = ['COMPRESS=LZW'])
#        for n in range(1,4):
#            outDataset.GetRasterBand(n).WriteArray(output_img[...,n-1])
#            outDataset.GetRasterBand(n).SetNoDataValue(0)
#
#    outDataset.SetProjection(source_proj)
#    outDataset.SetGeoTransform(source_trans)
#    del outDataset

def write_cog(output_file, data, source_trans, source_proj, metadata=None, nodata_value=-9999):
    """
    Writes a Cloud Optimized GeoTIFF (COG) file.

    Args:
        output_file (str): Path to the output COG file.
        data (numpy.ndarray): The spectral data to be written, with shape (rows, cols, bands).
        inds: (gdal dataset) GDAL dataset with metadata.
        metadata (dict, optional): Metadata to be added to the output file. Defaults to None.
        nodata_value (, optional): The nodata value to use. Defaults to -9999.
    """
    driver = gdal.GetDriverByName('MEM')

    od = data
    numpy_to_gdal = {
        np.dtype(np.float64): 7,
        np.dtype(np.float32): 6,
        np.dtype(np.int32): 5,
        np.dtype(np.uint32): 4,
        np.dtype(np.int16): 3,
        np.dtype(np.uint16): 2,
        np.dtype(np.uint8): 1,
    }

    ds = driver.Create('', od.shape[1], od.shape[0], od.shape[2], numpy_to_gdal[od.dtype])
    ds.SetGeoTransform(source_trans)
    ds.SetProjection(source_proj)
    if metadata:
        md = ds.GetMetadata()
        md.update(metadata)
        ds.SetMetadata(md)  

    for i in range(od.shape[2]):
        ds.GetRasterBand(i+1).WriteArray(od[:, :, i])
        ds.GetRasterBand(i+1).SetNoDataValue(nodata_value)
    ds.BuildOverviews('NEAREST', [2, 4, 8, 16, 32, 64, 128])
    
    tiff_driver = gdal.GetDriverByName('GTiff')
    output_dataset = tiff_driver.CreateCopy(output_file, ds, options=['COMPRESS=LZW', 'BIGTIFF=YES','COPY_SRC_OVERVIEWS=YES', 'TILED=YES', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256'])
    ds = None
    output_dataset = None

def plume_mask_threshold(input: np.array, pm, style='ch4'):

    y_locs = np.where(np.sum(pm > 0, axis=1))[0]
    x_locs = np.where(np.sum(pm > 0, axis=0))[0]

    plume_dat = input[y_locs[0]:y_locs[-1],x_locs[0]:x_locs[-1]].copy()
    plume_dat[pm[y_locs[0]:y_locs[-1],x_locs[0]:x_locs[-1]] == 0] = 0
    
    local_output_mask = None
    if style == 'ch4':
        plume_dat = gauss_blur(plume_dat, 3, preserve_nans=False)
        local_output_mask = plume_dat > 50
    else:
        plume_dat = gauss_blur(plume_dat, 5, preserve_nans=False)
        local_output_mask = plume_dat > 10000
    
    output_plume_mask = np.zeros(pm.shape,dtype=bool)
    output_plume_mask[y_locs[0]:y_locs[-1],x_locs[0]:x_locs[-1]] = local_output_mask
    
    return output_plume_mask



def spatial_temporal_filter(cov_df, coverage, roi, start_time, end_time):

    temporal_inds = np.where(np.logical_and(cov_df['properties.start_time'] >= pd.to_datetime(start_time) , 
                                            cov_df['properties.end_time'] <= pd.to_datetime(end_time) ))[0]
    spatial_inds = np.where(cov_df['geometry.coordinates'][:][temporal_inds].apply(lambda s,
                        roi=roi: s.intersects(roi)))[0]

    return [coverage['features'][i] for i in temporal_inds[spatial_inds]]

def roi_filter(coverage, roi):
    subdir = deepcopy(coverage) 
    cov_df = pd.json_normalize(coverage['features'])
    cov_df['geometry.coordinates'] = cov_df['geometry.coordinates'].apply(lambda s: Polygon(s[0]) )
    inds = np.where(cov_df['geometry.coordinates'].apply(lambda s,roi=roi: s.intersects(roi)))[0]
    subdir['features'] = [subdir['features'][i] for i in inds]
    return subdir
 
def time_filter(coverage, start_time, end_time):
    subdir = deepcopy(coverage)
    cov_df = pd.json_normalize(coverage['features'])
    inds = np.where(np.logical_and(pd.to_datetime(cov_df['properties.start_time']) >= pd.to_datetime(start_time) , pd.to_datetime(cov_df['properties.end_time']) <= pd.to_datetime(end_time)))[0]
    subdir['features'] = [subdir['features'][i] for i in inds]
    return subdir


class SerialEncoder(json.JSONEncoder):
    """Encoder for json to help ensure json objects can be passed to the workflow manager.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        else:
            return super(SerialEncoder, self).default(obj)


def add_fids(manual_annotations, coverage, manual_annotations_previous):
    manual_annotations_fid = deepcopy(manual_annotations)

    previous_plume_ids = []
    if manual_annotations_previous is not None:
        previous_plume_ids = [x['properties']['Plume ID'] for x in manual_annotations_previous['features']]

    updated_plumes=[]
    todel=[]

    # do some dataframe conversion once ahead of time to make things faster
    coverage_df = pd.json_normalize(coverage['features'])
    coverage_df['geometry.coordinates'] = coverage_df['geometry.coordinates'].apply(lambda s: Polygon(s[0]) )
    coverage_df['properties.start_time'] = pd.to_datetime(coverage_df['properties.start_time'])
    coverage_df['properties.end_time'] = pd.to_datetime(coverage_df['properties.end_time'])

    for _feat, feat in enumerate(manual_annotations['features']):
        logging.debug(f'Adding new fid {_feat} / {len(manual_annotations["features"])}')
        # If this key isn't present, then the full feature wasn't really added yet
        if 'R1 - Reviewed' not in feat['properties'].keys():
            todel.append(_feat)
            logging.info(f'R1 - Reviewed not in {feat["properties"]}')
            continue # This is insufficient
        plume_id = feat['properties']['Plume ID']

        if plume_id in previous_plume_ids:
            new_geom = feat['geometry']['coordinates']
            prev_idx = previous_plume_ids.index(plume_id)
            prev_geom = manual_annotations_previous['features'][prev_idx]['geometry']['coordinates']

            # check reviews
            rev_match = True
            for rl in ['R1 - Reviewed','R2 - Reviewed','R1 - VISIONS','R2 - VISIONS', 
                       'Psuedo-Origin', 'Sector', 'Sector Confidence', 'Time Range End', 'Time Range Start']:
                if rl in feat['properties'] and rl in manual_annotations_previous['features'][prev_idx]['properties']:
                    if feat['properties'][rl] != manual_annotations_previous['features'][prev_idx]['properties'][rl]:
                        rev_match = False

            if new_geom == prev_geom and rev_match:
                manual_annotations_fid['features'][_feat] = manual_annotations_previous['features'][prev_idx]
                #logging.debug(f'Geometries and properties the same in {feat["properties"]}...skipping safely')
                continue

        roi = Polygon(feat['geometry']['coordinates'][0])
        roi = shapely.buffer(roi, 0.01, join_style='mitre')
        subset_features = spatial_temporal_filter(coverage_df, coverage, roi, 
                                                  feat['properties']['Time Range Start'] + 'Z', 
                                                  feat['properties']['Time Range End'] + 'Z') 
        if len(subset_features) == 0:
            todel.append(_feat)
        else:
            fids = [subset_features[x]['properties']['fid'].split('_')[0] for x in range(len(subset_features))]
            manual_annotations_fid['features'][_feat]['properties']['fids'] = fids
            updated_plumes.append(_feat)

    if len(todel) > 0:
        logging.warning('Bad metadata for the following plumes:')
        for td in np.array(todel)[::-1]:
            msg = f'{manual_annotations_fid["features"][td]["properties"]["Plume ID"]}'
            logging.warning(msg)
        for td in np.array(todel)[::-1]:
            msg = f'Deleting entry due to bad metadata - check input {manual_annotations_fid["features"][td]["properties"]}'
            logging.warning(msg)
            manual_annotations_fid['features'].pop(td)

    updated_plumes = np.array([x for x in updated_plumes if x not in todel]) # shouldn't be necessary anymore, deosn't hurt
    for td in np.array(todel)[::-1]:
        updated_plumes[updated_plumes >= td] -= 1
    updated_plumes = updated_plumes.tolist()

    return manual_annotations_fid, updated_plumes


def add_orbits(annotations, indices_to_update, database):

    ind_to_pop = []
    update_ind_to_pop = []
    for _ind, ind in enumerate(indices_to_update):
        db_ret = [database.find_acquisition_by_id(fid) for fid in annotations['features'][ind]['properties']['fids']]
        orbits = [db_ret[_fid]['orbit'] for _fid, fid in enumerate(annotations['features'][ind]['properties']['fids'])]
        dcids  = [db_ret[_fid]['associated_dcid']  for _fid, fid in enumerate(annotations['features'][ind]['properties']['fids'])]

        scene_numbers  = [db_ret[_fid]['daac_scene'] if 'daac_scene' in db_ret[_fid].keys() else None for _fid, fid in enumerate(annotations['features'][ind]['properties']['fids']) ]

        un_orbits = np.unique(orbits)
        un_dcids = np.unique(dcids)

        if len(db_ret) == 0 or None in scene_numbers:
            logging.info(f'No FIDs or DAAC Scenes at {annotations["features"][ind]["properties"]["Plume ID"]}...skipping')
            annotations['features'][ind]['properties']['orbit'] = []
            ind_to_pop.append(ind)
            update_ind_to_pop.append(_ind)
            continue

        if len(un_dcids) > 1:
            logging.error(f'Ack - entry {annotations["features"][ind]} spans two dcids')
            annotations['features'][ind]['properties']['orbit'] = []
            ind_to_pop.append(ind)
            update_ind_to_pop.append(_ind)
            continue

        annotations['features'][ind]['properties']['orbit'] = un_orbits[0]
        annotations['features'][ind]['properties']['dcid'] = un_dcids[0]
        annotations['features'][ind]['properties']['daac_scenes'] = scene_numbers
    
    if len(ind_to_pop) > 0:
        logging.info('Bad plume list:')
    for ind in np.array(ind_to_pop)[::-1]:
        logging.info(annotations['features'][ind]['properties']['Plume ID'])
        annotations['features'].pop(ind)

    indices_to_update = np.array(indices_to_update)
    for _ind in np.array(update_ind_to_pop)[::-1]:
        indices_to_update[_ind:] -= 1
    indices_to_update = indices_to_update.tolist()

    for _ind in np.array(update_ind_to_pop)[::-1]:
        indices_to_update.pop(_ind)
    return annotations, indices_to_update

def write_color_plume(rawdat, plumes_mask, glt_ds, outname: str, style = 'ch4'):

    dat = rawdat.copy()
    colorized = np.zeros((dat.shape[0],dat.shape[1],3))

    dat[np.logical_not(plumes_mask)] = 0

    if style == 'ch4':
        dat /= 1500
        dat[dat > 1] = 1
        dat[dat < 0] = 0
        dat[dat == 0] = 0.01
        colorized[plumes_mask,:] = plt.cm.plasma(dat[plumes_mask])[...,:3]
    else:
        #dat /= 85000
        dat /= 85000
        dat[dat > 1] = 1
        dat[dat < 0] = 0
        dat[dat == 0] = 0.01
        colorized[plumes_mask,:] = plt.cm.viridis(dat[plumes_mask])[...,:3]


    colorized = (colorized * 255).astype(np.uint8)
    colorized[plumes_mask,:] = np.maximum(1, colorized[plumes_mask,:])

    write_cog(outname, colorized, glt_ds, nodata_value=0)


def write_color_quicklook(indat, output_file):
    dat = indat.copy()
    mask = dat != -9999
    dat[dat < 0] = 0
    dat = dat /1500.
    output = np.zeros((indat.shape[0],indat.shape[1],3),dtype=np.uint8)
    output[mask,:] = np.round(plt.cm.plasma(dat[mask])[...,:3] * 255).astype(np.uint8)
    output[mask,:] = np.maximum(1, output[mask])


    memdriver = gdal.GetDriverByName('MEM')
    memdriver.Register()
    outDataset = memdriver.Create('',dat.shape[1],dat.shape[0],3,gdal.GDT_Byte)
    for n in range(1,4):
        outDataset.GetRasterBand(n).WriteArray(output[...,n-1])
        outDataset.GetRasterBand(n).SetNoDataValue(0)

    driver = gdal.GetDriverByName('PNG')
    driver.Register()
    dst_ds = driver.CreateCopy(output_file, outDataset, strict=0)
    del dst_ds, outDataset
 

def write_json(output_file, plume_dict, scene_names):
    outdict = {"crs": {"properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84" }, "type": "name"},"features":[],"name":"methane_metadata","type":"FeatureCollection" }
    outdict['features'].append(deepcopy(plume_dict))

    outdict['features'][0]['properties']['DAAC Scene Names'] = scene_names
    del outdict['features'][0]['properties']['style']
    del outdict['features'][0]['properties']['Data Download']

    with open(output_file, 'w') as fout:
        fout.write(json.dumps(outdict, cls=SerialEncoder)) 



def prep_predictor_image(predictor, data, ptype):
    image = data.copy()
    if ptype == 'co2':
        image = (image / 1500).astype(np.uint8)
    else:
        image = (image / 100000).astype(np.uint8)
    image[image > 1] = 1
    image[image < 0] = 0
    image *= 255
    oi = np.zeros((image.shape[0],image.shape[1],3),dtype=np.uint8)
    oi[...] = image.squeeze()[:,:,np.newaxis]

    predictor.set_image(oi)

def rawspace_coordinate_conversion(glt, coordinates, trans, ortho=False):
    rawspace_coords = []
    for ind in coordinates:
        glt_ypx = int(round((ind[1] - trans[3])/ trans[5]))
        glt_xpx = int(round((ind[0] - trans[0])/ trans[1]))
        if ortho:
            rawspace_coords.append([glt_xpx,glt_ypx])
        else:
            lglt = glt[glt_ypx, glt_xpx,:]
            rawspace_coords.append(lglt.tolist())
    return rawspace_coords


def sam_segmentation(predictor, rawspace_coords, manual_mask, n_input_points=20):
    min_x = np.min([x[0] for x in rawspace_coords])
    max_x = np.max([x[0] for x in rawspace_coords])
    min_y = np.min([x[1] for x in rawspace_coords])
    max_y = np.max([x[1] for x in rawspace_coords])
    bbox = np.array([min_x, min_y, max_x, max_y])

    input_labels, input_points = [],[]
    for _n in range(n_input_points):
        pt = [np.random.randint(min_x, max_x), np.random.randint(min_y,max_y)]
        input_labels.append(manual_mask[pt[1],pt[0]])
        input_points.append(pt)
        
    masks, _, _ = predictor.predict(
                                    point_coords=np.array(input_points),
                                    point_labels=np.array(input_labels),
                                    box=bbox,
                                    multimask_output=False,
                                   )

    return masks[0,...]

def trim_plume(plume_mask, trans, nodata_value=0):
    y_locs = np.where(np.sum(plume_mask != nodata_value, axis=1))[0]
    x_locs = np.where(np.sum(plume_mask != nodata_value, axis=0))[0]

    if len(y_locs) == 0 or len(x_locs) == 0:
        raise ValueError('No valid plume pixels found in mask to trim')

    trimmed_mask = plume_mask[y_locs[0]:y_locs[-1]+1, x_locs[0]:x_locs[-1]+1].copy()
    outtrans = list(deepcopy(trans))
    outtrans[0] += x_locs[0] * trans[1]
    outtrans[3] += y_locs[0] * trans[5]
    return trimmed_mask, outtrans

def build_plume_properties(plume_input_properties, plume_geometry, plume_data, trans, data_version, xsize_m, ysize_m, nodata_value=-9999):

    loc_pp = deepcopy(plume_input_properties)
    loc_pp['geometry'] = plume_geometry

    props = {}
    props['Plume ID'] = loc_pp['Plume ID']
    props['Scene FIDs'] = loc_pp['fids']
    props['Orbit'] = loc_pp['orbit']
    props['DCID'] = loc_pp['dcid']
    props['DAAC Scene Numbers'] = loc_pp['daac_scenes']
    props["Data Download"] = get_daac_link(props, data_version)


    inplume_dat = plume_data[plume_data != nodata_value]

    start_datetime = datetime.datetime.strptime(loc_pp['fids'][0][4:], "%Y%m%dt%H%M%S")
    end_datetime = start_datetime + datetime.timedelta(seconds=1)

    start_datetime = start_datetime.strftime("%Y-%m-%dT%H:%M:%SZ")
    end_datetime = end_datetime.strftime("%Y-%m-%dT%H:%M:%SZ")

   
    
    maxval = np.nanmax(inplume_dat)
    if len(inplume_dat) < 5 or maxval < 200 or np.isnan(maxval):
        logging.warning(f'Plume {loc_pp["Plume ID"]} has insufficient data...skipping')
        return None

    rawloc = np.where(plume_data == maxval)
    maxval = np.round(maxval)

    #sum and convert to kg.  conversion:
    # ppm m / 1e6 ppm * x_pixel_size(m)*y_pixel_size(m) 1e3 L / m^3 * 1 mole / 22.4 L * 0.01604 kg / mole
    ime_scaler = (1.0/1e6)* ((np.abs(xsize_m*ysize_m))/1.0) * (1000.0/1.0) * (1.0/22.4)*(0.01604/1.0)

    ime = np.nansum(inplume_dat) * ime_scaler
    ime_p = np.nansum(inplume_dat[inplume_dat > 0]) * ime_scaler

    ime = np.round(ime,2)
    #ime_uncert = np.round(len(inplume_dat) * background * ime_scaler,2)

    max_loc_y = trans[3] + trans[5]*(rawloc[0][0]+0.5)
    max_loc_x = trans[0] + trans[1]*(rawloc[1][0]+0.5)

    props["UTC Time Observed"] = start_datetime
    props["map_endtime"] = end_datetime
    props["Max Plume Concentration (ppm m)"] = maxval
    #props["Concentration Uncertainty (ppm m)"] = background
    #props["Integrated Methane Enhancement (kg CH4)"] = ime
    #props["Integrated Methane Enhancement - Positive (kg CH4)"] = ime_p
    #props["Integrated Methane Enhancement Uncertainty (kg CH4)"] = ime_uncert
    props["Latitude of max concentration"] = np.round(max_loc_y, 5)
    props["Longitude of max concentration"] = np.round(max_loc_x, 5)
    

    # For R1 Review
    if not loc_pp['R1 - Reviewed']:
        props['style'] = {"maxZoom": 20, "minZoom": 0, "color": "red", 'opacity': 1, 'weight': 2, 'fillOpacity': 0}
    
    # For R2 Review
    if loc_pp['R1 - Reviewed'] and loc_pp['R1 - VISIONS'] and not loc_pp['R2 - Reviewed']:
        props['style'] = {"maxZoom": 20, "minZoom": 0, "color": "green", 'opacity': 1, 'weight': 2, 'fillOpacity': 0}

    # Accept
    if loc_pp['R1 - Reviewed'] and loc_pp['R1 - VISIONS'] and loc_pp['R2 - Reviewed'] and loc_pp['R2 - VISIONS']:
        props['style'] = {"maxZoom": 20, "minZoom": 0, "color": "white", 'opacity': 1, 'weight': 2, 'fillOpacity': 0}
    
    # Reject
    if (loc_pp['R1 - Reviewed'] and not loc_pp['R1 - VISIONS']) or (loc_pp['R2 - Reviewed'] and not loc_pp['R2 - VISIONS']):
        props['style'] = {"maxZoom": 20, "minZoom": 0, "color": "yellow", 'opacity': 1, 'weight': 2, 'fillOpacity': 0}

    poly_res = {"geometry": loc_pp['geometry'],
               "type": "Feature",
               "properties": props}

    props['style']['radius'] = 10
    point_res = {"geometry": {"coordinates": [max_loc_x, max_loc_y, 0.0], "type": "Point"},
               "properties": props,
               "type": "Feature"}
    
    return poly_res, point_res


def get_daac_link(feature, product_version):
    prod_v = product_version.split('V')[-1]
    fid=feature['Scene FIDs'][0]
    cid= feature['Plume ID'].split('-')[-1].zfill(6)
    link = f'https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2BCH4PLM.{prod_v}/EMIT_L2B_CH4PLM_{prod_v}_{fid[4:12]}T{fid[13:19]}_{cid}/EMIT_L2B_CH4PLM_{prod_v}_{fid[4:12]}T{fid[13:19]}_{cid}.tif'
    return link


def get_sds_cog(fid, enh_version, dtype='ch4'):
    date = fid[4:12]
    path=f'/store/emit/ops/data/acquisitions/{date}/{fid.split("_")[0]}/ghg/{dtype}/{fid.split("_")[0]}*_ghg_ort{dtype}_b0106_{enh_version}.tif'
    return path


def write_geojson_linebyline(output_file, outdict):
    with open (output_file, 'w') as f:
        # Write opening brace
        f.write('{\n')
        
        # Write all properties except 'features' with 4-space indent
        properties = {k: v for k, v in outdict.items() if k != 'features'}
        for i, (key, value) in enumerate(properties.items()):
            value_str = json.dumps(value, indent=4, cls=SerialEncoder)
            # Indent the value if it's multi-line
            value_lines = value_str.split('\n')
            if len(value_lines) > 1:
                indented_value = '\n'.join(['    ' + line if idx > 0 else line for idx, line in enumerate(value_lines)])
            else:
                indented_value = value_str
            f.write(f'    "{key}": {indented_value},\n')
        
        # Write features array opening
        f.write('    "features": [\n')
        
        # Write each feature line by line
        for i, feat in enumerate(outdict['features']):
            line = json.dumps(feat, separators=(",", ":"), cls=SerialEncoder)
            if i < len(outdict['features']) - 1:
                f.write(line + ",\n")
            else:
                f.write(line + "\n")
        
        # Close features array and main object
        f.write('    ]\n')
        f.write('}')

def global_metadata(data_version, software_version=None):
    extra_metadata = {}
    if software_version:
        extra_metadata['software_build_version'] = software_version
    else:
        cmd = ["git", "symbolic-ref", "-q", "--short", "HEAD", "||", "git", "describe", "--tags", "--exact-match"]
        output = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
        if output.returncode != 0:
            raise RuntimeError(output.stderr.decode("utf-8"))
        extra_metadata['software_build_version'] = output.stdout.decode("utf-8").replace("\n", "")

    extra_metadata['product_version'] = data_version
    extra_metadata['keywords'] = "Imaging Spectroscopy, minerals, EMIT, dust, radiative forcing"
    extra_metadata['sensor'] = "EMIT (Earth Surface Mineral Dust Source Investigation)"
    extra_metadata['instrument'] = "EMIT"
    extra_metadata['platform'] = "ISS"
    extra_metadata['Conventions'] = "CF-1.63"
    extra_metadata['institution'] = "NASA Jet Propulsion Laboratory/California Institute of Technology"
    extra_metadata['license'] = "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
    extra_metadata['naming_authority'] = "LPDAAC"
    extra_metadata['date_created'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    extra_metadata['keywords_vocabulary'] = "NASA Global Change Master Directory (GCMD) Science Keywords"
    extra_metadata['stdname_vocabulary'] = "NetCDF Climate and Forecast (CF) Metadata Convention"
    extra_metadata['creator_name'] = "Jet Propulsion Laboratory/California Institute of Technology"
    extra_metadata['creator_url'] = "https://earth.jpl.nasa.gov/emit/"
    extra_metadata['project'] = "Earth Surface Mineral Dust Source Investigation"
    extra_metadata['project_url'] = "https://earth.jpl.nasa.gov/emit/"
    extra_metadata['publisher_name'] = "NASA LPDAAC"
    extra_metadata['publisher_url'] = "https://lpdaac.usgs.gov"
    extra_metadata['publisher_email'] = "lpdaac@usgs.gov"
    extra_metadata['identifier_product_doi_authority'] = "https://doi.org"
    extra_metadata['title'] = "EMIT"
    extra_metadata['Units']= 'ppm m'
    return extra_metadata


def get_metadata(plume_dict, global_metadata):
    scene_names = []
    for _s in range(len(plume_dict['properties']['Scene FIDs'])):
        fid =plume_dict['properties']['Scene FIDs'][_s]
        scene =plume_dict['properties']['DAAC Scene Numbers'][_s]
        orbit =plume_dict['properties']['Orbit']
        scene_names.append(f'EMIT_L2B_CH4ENH_{global_metadata["product_version"]}_{fid[4:12]}T{fid[13:19]}_{orbit}_{scene}')

    metadata = {
            'Plume_Complex': plume_dict['properties']['Plume ID'],
            #'Estimated_Uncertainty_ppmm': plume_dict['properties']['Concentration Uncertainty (ppm m)'],
            'UTC_Time_Observed': plume_dict['properties']['UTC Time Observed'],
            #Source_Scenes - match full conventions 
            'Source_Scenes': ','.join(scene_names),
            'Latitude of max concentration': plume_dict['properties']['Latitude of max concentration'],
            'Longitude of max concentration': plume_dict['properties']['Longitude of max concentration'],
            'Max Plume Concentration (ppm m)': plume_dict['properties']['Max Plume Concentration (ppm m)'],
            }
    metadata.update(global_metadata)
    return metadata

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

    # Global File names
    annotation_file_raw = os.path.join(args.out_dir, "manual_annotation_raw.json") # Straight from MMGIS
    annotation_file = os.path.join(args.out_dir, "manual_annotation.json") # MMGIS + FIDs + Orbits 
    previous_annotation_file = os.path.join(args.out_dir, "previous_manual_annotation.json") # Last editted version of annotation file
    output_json = os.path.join(args.out_dir, 'combined_plume_metadata.json') # Output (public facing) metadata file
    delivery_dir = os.path.join(args.out_dir, 'delivery') # Delivery file directory
    os.makedirs(delivery_dir, exist_ok=True)


    for run in range(max_runs):
        logging.debug('Loading Data')
        ######## Step 1 ###########

        # Update annotations file load
        print_and_call(f'curl "https://popo.jpl.nasa.gov/mmgis/API/files/getfile" -H "Authorization:Bearer {args.key}" --data-raw "id={args.id}" > {annotation_file_raw}')
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
        manual_annotations, new_plumes = add_fids(manual_annotations, coverage, manual_annotations_previous)
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
        manual_annotations, new_plumes     = add_orbits(manual_annotations, new_plumes, database)
        logging.info(f'New plume, post orbit: {len(new_plumes)}')



        ######## Step 2 ###########
        unique_fids = np.unique([sublist for x in new_plumes for sublist in manual_annotations['features'][x]['properties']['fids']])
        unique_orbits = np.unique([manual_annotations['features'][x]['properties']['orbit'] for x in new_plumes]).tolist()
        unique_dcids = np.unique([manual_annotations['features'][x]['properties']['dcid'] for x in new_plumes]).tolist()

        # Dump out the udpated manual annotations set, so it holds FIDs / orbits for next round
        logging.debug('Dump out the updated manual annotations set, overwrting original')
        with open(annotation_file, 'w') as fout:
            fout.write(json.dumps(manual_annotations, cls=SerialEncoder)) 

        for feat in manual_annotations['features']:
            if 'dcid' not in feat['properties'].keys():
                print('nodcid:  ', feat)


        # Now step through each DCID, get the relevant FIDs, mosaic, and cut out each plume
        for dcid in unique_dcids[:1]:

            logging.info(f'Processing DCID {dcid}...')
            #plumes_idx_in_dcid = [x for x in range(len(manual_annotations['features'])) if manual_annotations['features'][x]['properties']['dcid'] == dcid]
            plumes_idx_in_dcid = [x for x in new_plumes if manual_annotations['features'][x]['properties']['dcid'] == dcid]
            fids_in_dcid = np.unique([sublist for x in plumes_idx_in_dcid for sublist in manual_annotations['features'][x]['properties']['fids']])
            logging.info(f'...found {len(plumes_idx_in_dcid)} plumes to update and {len(fids_in_dcid)} FIDs')

            ort_dat_files = [get_sds_cog(fid, args.enh_data_version, dtype=args.type) for fid in fids_in_dcid]
            dcid_ort_vrt_file = os.path.join(args.out_dir, f'dcid_{dcid}_mf_ort.vrt')
            dcid_ort_tif_file = os.path.join(args.out_dir, f'dcid_{dcid}_mf_ort.tif')
            print_and_call(f'gdalbuildvrt {dcid_ort_vrt_file} {" ".join(ort_dat_files)} --quiet')
            print_and_call(f'gdal_translate {dcid_ort_vrt_file} {dcid_ort_tif_file} -co COMPRESS=LZW --quiet')

            ort_ds = gdal.Open(dcid_ort_tif_file)
            ortdat = ort_ds.ReadAsArray().squeeze()
            trans = ort_ds.GetGeoTransform()

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

                rawspace_coords = rawspace_coordinate_conversion([], feat['geometry']['coordinates'][0], trans, ortho=True)
                manual_mask = rasterize(shapes=[Polygon(rawspace_coords)], out_shape=(ortdat.shape[0],ortdat.shape[1])) # numpy binary mask for manual IDs

                plumestyle = 'classic'
                if 'Delineation Mode' in feat['properties'].keys():
                    plumestyle = feat['properties']['Delineation Mode']

                loc_fid_mask = None # change name to dcid_plume_mask
                if plumestyle == 'classic':
                    loc_fid_mask = plume_mask_threshold(ortdat.copy(), manual_mask, style=args.type)
                elif plumestyle == 'manual':    
                    loc_fid_mask = manual_mask.astype(bool)
                
                if combined_mask is None:
                    combined_mask = loc_fid_mask.copy()
                else:
                    combined_mask = np.logical_or(combined_mask, loc_fid_mask)

                cut_plume_mask, newp_trans = trim_plume(loc_fid_mask, trans)
                tocut = ortdat.copy()
                tocut[loc_fid_mask == 0] = -9999
                cut_plume_data, _ = trim_plume(tocut, trans, nodata_value=-9999)



                ############  Step 4 ###########
                outbase = plume_working_basename(args.out_dir, feat)
                outmask_finepoly_file = os.path.join(outbase + '_finepolygon.json')
                outmask_poly_file = os.path.join(outbase + '_polygon.json')
                outmask_ort_file = os.path.join(outbase + '_mask_ort.tif')

                # Write mask file and save tif reference
                #write_output_file(newp_trans, ort_ds.GetProjection(), cut_plume_mask, outmask_ort_file)
                #write_output_file(newp_trans, ort_ds.GetProjection(), cut_plume_data, outmask_ort_file)
                write_cog(outmask_ort_file, cut_plume_mask.reshape((cut_plume_mask.shape[0], cut_plume_mask.shape[1],1)).astype(np.uint8), 
                          newp_trans, ort_ds.GetProjection(), nodata_value=0)
                
                if os.path.isfile(outmask_poly_file):
                    os.remove(outmask_poly_file)
                if os.path.isfile(outmask_finepoly_file):
                    os.remove(outmask_finepoly_file)
                print_and_call(f'gdal_polygonize {outmask_ort_file} {outmask_finepoly_file} -f GeoJSON -mask {outmask_ort_file} -8 -quiet')
                print_and_call(f'ogr2ogr {outmask_poly_file} {outmask_finepoly_file} -f GeoJSON -lco RFC7946=YES -simplify {trans[1]} --quiet')

                ## Save a running tile of the ort files and polygon files for this DCID, only for mosaicing
                #dcid_mask_tif_files.append(outmask_ort_file)
                #dcid_mask_poly_files.append(outmask_poly_file)

                # Read back in the polygon we just wrote
                plume_to_add = json.load(open(outmask_poly_file))['features']
                if len(plume_to_add) > 1:
                    logging.warning(f'ACK - multiple polygons from one Plume ID in file {outmask_poly_file}')
                plume_to_add[0]['geometry']['coordinates'] = [[np.round(x[0],5), np.round(x[1],5)] for x in plume_to_add[0]['geometry']['coordinates'][0]]

                poly_plume, point_plume = build_plume_properties(feat['properties'], plume_to_add[0]['geometry'], cut_plume_data, 
                                                                 newp_trans, args.data_version, xsize_m, ysize_m)
                updated_plumes_point.append(point_plume)
                updated_plumes_poly.append(poly_plume)


                # These are the real delivery files
                delivery_base = plume_delivery_basename(delivery_dir, poly_plume)
                os.makedirs(os.path.dirname(delivery_base), exist_ok=True)

                delivery_raster_file = os.path.join(delivery_base + '.tif')
                delivery_ql_file = os.path.join(delivery_base + '.png')
                delivery_json_file = os.path.join(delivery_base + '.json')

                meta = get_metadata(poly_plume, global_metadata(data_version=args.data_version))
                write_cog(delivery_raster_file, cut_plume_data.reshape((cut_plume_data.shape[0], cut_plume_data.shape[1],1)).astype(np.uint8), 
                          newp_trans, ort_ds.GetProjection(), nodata_value=-9999, metadata=meta)
                write_color_quicklook(cut_plume_data, delivery_ql_file)
                write_json(delivery_json_file, poly_plume, meta['Source_Scenes'])

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

        write_geojson_linebyline(output_json, outdict) # Final write


        # Sync
        print_and_call(f'cp {annotation_file} {previous_annotation_file}')
        #print_and_call(f'rsync -a --info=progress2 {tile_dir}/{od_date}/ brodrick@$EMIT_SCIENCE_IP:/data/emit/mmgis/mosaics/{args.type}_plume_tiles_working/{od_date}/ --delete')
        #print_and_call(f'rsync {output_json} brodrick@$EMIT_SCIENCE_IP:/data/emit/mmgis/coverage/converted_manual_{args.type}_plumes.json')





        # Estimating the background would require going through *ALL* plumes in the DCID, not just the new ones - no longer needed, just leaving as a reminder
        #background = ortdat[np.logical_and(combined_mask == 0, ortdat != -9999)]















 
if __name__ == '__main__':
    main()



