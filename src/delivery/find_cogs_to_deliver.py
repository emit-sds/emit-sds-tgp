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
# Author: Winston Olson-Duvall, winston.olson-duvall@jpl.nasa.gov

import argparse
import glob
import json
import os
import requests


def main():
    # Set up args
    parser = argparse.ArgumentParser(description="Deliver GHG products to LP DAAC")
    # parser.add_argument("ghg", help="Must specify either 'ch4' or 'co2'")
    parser.add_argument("base_dir", help="Base directory - use this to determine which GHG")
    parser.add_argument("--cmr", action="store_true", help="Use CMR query to identify already published granules")
    parser.add_argument("--token_file", help="Path to text file containing CMR token")
    args = parser.parse_args()

    # Default base_dirs in the past have been
    # - /scratch/brodrick/methane/visions_delivery_20240409
    # - /scratch/brodrick/methane/visions_delivery_co2
    # - /store/brodrick/repos/emit-sds-tgp/V002_ch4/daac/

    base_dir = args.base_dir
    ghg = "ch4"
    if "co2" in base_dir:
        ghg = "co2"

    if not args.cmr:
        # If not querying CMR, then check the filesystem for cmr.json files
        cogs = glob.glob(os.path.join(base_dir, "20*/*tif"))
        cmr_files = glob.glob(os.path.join(base_dir, "20*/*cmr.json"))
        basenames = [os.path.basename(f).replace(".cmr.json", "") for f in cmr_files]
        for cog in cogs:
            if os.path.basename(cog).replace(".tif", "") not in basenames:
                print(cog)
    else:
        # If querying CMR, then need compare COGs against CMR records
        CMR_OPS = 'https://cmr.earthdata.nasa.gov/search'  # CMR API Endpoint
        url = f'{CMR_OPS}/{"granules"}'

        # C2748088093-LPCLOUD CH4PLM emit20220820t101039_CH4_PlumeComplex-2715.tif EMIT_L2B_CH4PLM_001_20220820T101039_002715
        # C2867824144-LPCLOUD CO2PLM emit20231204t082834_CO2_PlumeComplex-277.tif EMIT_L2B_CO2PLM_001_20231204T082834_000277

        # V001 - kept here for reference 
        # plm_coll = "C2748088093-LPCLOUD" if ghg == "ch4" else "C2867824144-LPCLOUD"

        # V002
        plm_coll = "C3242707413-LPCLOUD" if ghg == "ch4" else "C3244277342-LPCLOUD"
        datetime_range = '2022-08-01T00:00:00Z,2027-01-01T00:00:00Z'  # Overall date range of granules to be searched
        page_size = 2000

        # Get token from token file
        if args.token_file is not None:
            token_file = args.token_file
        else:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            token_file = os.path.join(script_dir, "cmr_token.txt")
        with open(token_file, "r") as f:
            token = f.read().replace("\n", "")

        # Check the plume collection
        # print(f"Checking {ghg} plume collection {plm_coll}")
        response = requests.get(url,
                                params={'concept_id': plm_coll,
                                        'temporal': datetime_range,
                                        'page_size': page_size
                                        },
                                headers={
                                    'Accept': 'application/json',
                                    'Authorization': f'Bearer {token}'
                                }
                                )
        # print(response.status_code)
        # print(f"Number of granules found: {response.headers['CMR-Hits']}")  # Resulting quantity of granules/items.
        granules = response.json()['feed']['entry']
        results = [g["title"].split("_")[4] + "_" + g["title"].split("_")[5] for g in granules]
        plm_cogs = glob.glob(os.path.join(base_dir, "20*/*tif"))
        # print(f"Number of plm COGs on filesystem: {len(plm_cogs)}")
        for cog in plm_cogs:
            # If cog not in cmr results, then print (use timestamp_plumeid to check unique)
            timestamp = os.path.basename(cog)[4:19].upper()
            plume_id = os.path.basename(cog).split("-")[-1].replace(".tif", "").zfill(6)
            unique_plm = f"{timestamp}_{plume_id}"
            if unique_plm not in results:
                print(cog)


if __name__ == '__main__':
    main()