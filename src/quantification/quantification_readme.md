

# Quantifying Emissions
Computing the emission rate for a given plume requires estimating the Integrated Mass Enhancement (IME), the plume length (fetch), and the windspeed. The first two are computed here using compute_flux.py while the third is found using HRRR_and_ERA5_windspeeds.py. These two codes produce csv outputs that are combined using compute_Q_and_Qunc_simple_IME.py to generate a final csv output file that contains the estimated emission rate and related uncertainties.

## Required Inputs
compute_flux.py requires both the concentration length estimates and corresponding sesitivity and uncertainty (https://github.com/emit-sds/emit-ghg and downloadable from https://www.earthdata.nasa.gov/data/instruments/emit-imaging-spectrometer). It also requires a "pseudo-origin" marking the approximate location of the plume origin and plume boundary polygon.

## IME and Fetch
compute_flux.py computes the IME and plume fetch. An example of how to call compute_flux.py is found in simpleIME.py:

```
python src/quantification/simpleIME.py --data_folder /store/emit/ops/data/acquisitions/ --output_path OUTPUT_PATH --instrument emit --plume_json_filename /store/brodrick/methane/ch4_plumedir/previous_manual_annotation_oneback.json
```
You can also specify a subset of plumes to run using --plume_id id1 id2 ...
For example: 

```
python src/quantification/simpleIME.py --data_folder /store/emit/ops/data/acquisitions/ --output_path OUTPUT_PATH --instrument emit --plume_json_filename /store/brodrick/methane/ch4_plumedir/previous_manual_annotation_oneback.json --plume_id CH4_PlumeComplex-10 CH4_PlumeComplex-11
```

## Windspeed

Windspeed estimates are collected from HRRR and ERA-5 reanalysis data using HRRR_and_ERA5_windspeeds.py. An example call is:

```
python src/quantification/windspeeds.py /store/jfahlen/EMIT_met_data/EMIT_wind_data_HRRR_ERA5_fromCDS_0000.csv --write-rate 1 --plume-file /store/brodrick/methane/ch4_plumedir/previous_manual_annotation_oneback.json
```

Note that the first parameter is an input file. The output of the code will have the same name but with the number incremented. If the file does not exist, then it will be created. If it does exist, then new outputs will be appended. For the example above, if the file in the first parameter does not exist, then the code will store the windspeeds for all the plumes in --plume-file in EMIT_wind_data_HRRR_ERA5_fromCDS_0001.csv (note the increment). If it did exist, say from a previous run of the program, then EMIT_wind_data_HRRR_ERA5_fromCDS_0001.csv will contain all the plumes in the original plus the windspeeds from the additional plumes listed in --plume-file.

# Combining Outputs To Yield Emission Rate

The outputs of the two steps above are combined using compute_Q_and_Qunc_simple_IME.py. An example call is

```
python compute_Q_and_Qunc_simple_IME.py OUTPUT_from_Simple_IME.csv OUTPUT_from_windspeed.csv emission_rates.csv
```

See the --help outputs for details.