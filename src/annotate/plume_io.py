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

import numpy as np
import json
from osgeo import gdal
from copy import deepcopy
import subprocess
import matplotlib.pyplot as plt
import datetime
from quantification import compute_Q_and_uncertainty_utils


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
 

def write_delivery_json(output_file, plume_dict, scene_names, deliver_emissions, indent=4):
    outdict = {"crs": {"properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84" }, "type": "name"},"features":[],"name":"methane_metadata","type":"FeatureCollection" }
    outdict['features'].append(deepcopy(plume_dict))
    authorized_keys = ['Plume ID', 'Orbit', 'DCID', 'DAAC Scene Names', 'UTC Time Observed', 
                       'Max Plume Concentration (ppm m)', 'Latitude of max concentration', 'Longitude of max concentration', 'Wind Speed (m/s)', 
                       'Wind Speed Std (m/s)', 'Wind Speed Source', 'Emissions Rate Estimate (kg/hr)', 'Emissions Rate Estimate Uncertainty (kg/hr)', 'Fetch Length (m)']

    for key in list(outdict['features'][0]['properties'].keys()):
        if key not in authorized_keys:
            del outdict['features'][0]['properties'][key]
    if not deliver_emissions:
        outdict['features'][0]['properties'].update(compute_Q_and_uncertainty_utils.EMISSIONS_DELIVERY_INFO)

    outdict['features'][0]['properties']['DAAC Scene Names'] = scene_names

    with open(output_file, 'w') as fout:
        fout.write(json.dumps(outdict, indent=indent, cls=SerialEncoder)) 


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

def write_external_geojson(input_file, output_file):
    valid_keys = ['Plume ID', 'Scene FIDs', 'Orbit', 'DCID', 'DAAC Scene Names', 'Data Download', 'UTC Time Observed', 
                  'Max Plume Concentration (ppm m)', 'Latitude of max concentration', 'Longitude of max concentration', 'style', 
                  'Wind Speed (m/s)', 'Wind Speed Std (m/s)', 'Wind Speed Source', 'Emissions Rate Estimate (kg/hr)', 
                  'Emissions Rate Estimate Uncertainty (kg/hr)', 'Fetch Length (m)']
    outdict = json.load(open(input_file, 'r'))
    for _f in range(len(outdict['features'])):
        for key in list(outdict['features'][_f]['properties'].keys()):
            if key not in valid_keys:
                del outdict['features'][_f]['properties'][key]
    write_geojson_linebyline(output_file, outdict)



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
        scene =plume_dict['properties']['daac_scenes'][_s]
        orbit =plume_dict['properties']['Orbit']
        scene_names.append(f'EMIT_L2B_CH4ENH_{global_metadata["product_version"]}_{fid[4:12]}T{fid[13:19]}_{orbit}_{scene}')

    metadata = {
            'Plume_Complex': plume_dict['properties']['Plume ID'],
            #'Estimated_Uncertainty_ppmm': plume_dict['properties']['Concentration Uncertainty (ppm m)'],
            'UTC_Time_Observed': plume_dict['properties']['UTC Time Observed'],
            #Source_Scenes - match full conventions 
            'DAAC Scene Names': scene_names,
            'Latitude of max concentration': plume_dict['properties']['Latitude of max concentration'],
            'Longitude of max concentration': plume_dict['properties']['Longitude of max concentration'],
            'Max Plume Concentration (ppm m)': plume_dict['properties']['Max Plume Concentration (ppm m)'],
            }
    metadata.update(global_metadata)
    return metadata


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


def trim_plume(plume_mask_in, trans, badmask=None, nodata_value=0, buffer=0):

    plume_mask = plume_mask_in
    if badmask is not None:
        plume_mask = plume_mask_in.copy()
        plume_mask[badmask] = nodata_value

    y_locs = np.where(np.sum(plume_mask != nodata_value, axis=1))[0]
    x_locs = np.where(np.sum(plume_mask != nodata_value, axis=0))[0]

    if len(y_locs) == 0 or len(x_locs) == 0:
        return None, None
        #raise ValueError('No valid plume pixels found in mask to trim')
    
    min_y = max(y_locs[0] - buffer, 0)
    max_y = min(y_locs[-1] + buffer, plume_mask.shape[0] - 1)
    min_x = max(x_locs[0] - buffer, 0)
    max_x = min(x_locs[-1] + buffer, plume_mask.shape[1] - 1)

    trimmed_mask = plume_mask[min_y:max_y+1, min_x:max_x+1].copy()
    outtrans = list(deepcopy(trans))
    outtrans[0] += min_x * trans[1]
    outtrans[3] += min_y * trans[5]
    return trimmed_mask, outtrans

def get_window(rawspace_coords, trans, ds_size, buffer_px):

    xs = [x[0] for x in rawspace_coords]
    ys = [x[1] for x in rawspace_coords]
    min_x = max(0, int(np.min(xs)) - buffer_px)
    max_x = min(ds_size[1] - 1, int(np.max(xs)) + buffer_px)
    min_y = max(0, int(np.min(ys)) - buffer_px)
    max_y = min(ds_size[0] - 1, int(np.max(ys)) + buffer_px)
    
    win_w = max_x - min_x
    win_h = max_y - min_y
    
    if win_w <= 0 or win_h <= 0:
        return None, None, None
    
    # Window transform
    window_trans = list(trans)
    window_trans[0] += min_x * trans[1]
    window_trans[3] += min_y * trans[5]

    # Local coords for mask
    local_coords = [[x[0] - min_x, x[1] - min_y] for x in rawspace_coords]

    return (min_x, min_y, win_w, win_h), window_trans, local_coords
