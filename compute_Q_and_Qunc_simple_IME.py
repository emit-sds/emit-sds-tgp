#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import click
import compute_Q_and_uncertainty_utils


@click.command()
@click.argument('simple_ime_csv_filename', type=click.Path(file_okay=True), required=True)
@click.argument('wind_csv_filename',    type=click.Path(file_okay=True), required=True)
@click.argument('output_csv_filename', type=click.Path(file_okay=True), required=True)

# Options that a user may modify
@click.option('--fetch_unc_frac', type=float, default=0.,  show_default=True, help='Fetch length uncertainty fraction')
@click.option('--c_q_key',        default='C_Q_mask_kg_hr_mpers',  show_default=True, help='C_Q column key')
@click.option('--fetch_key',      default='fetch',                 show_default=True, help='fetch column key')
@click.option('--sum_c2_key',     default='sum_C_sqrd',            show_default=True, help='Squared pixel mass column key')

# Don't change unless wind speed keys have been added
@click.option('--primary_windspeed_source',   default='hrrr', show_default=True, help='Primary wind-speed source name')
@click.option('--secondary_windspeed_source', default='era5', show_default=True, help='Secondary wind-speed source name')

def compute_Q_and_Qunc_simple_IME(
        simple_ime_csv_filename, wind_csv_filename, output_csv_filename,
        fetch_unc_frac, c_q_key, fetch_key, sum_c2_key,
        primary_windspeed_source, secondary_windspeed_source):
    """
    Combine Simple-IME emissions with wind data to compute Q (kg/hr) and its
    uncertainty.  Only the fetch-uncertainty fraction is intended to be changed
    by most users; wind-speed columns are generated automatically from the
    two source names supplied above. The secondary windspeed source will be used
    only when the primary is not available.
    
        The intent is to take the C_Q values estimated by Simple IME and combine them with the windspeeds to calculate:
        - the emission rate Q in kg / hr
        - the emission rate uncertainty Q_unc in kg / hr
        - the components of Q_unc

    The fetch length uncertainty is a user-defined quantity provided as a fraction of the fetch length.
    For example, if fetch_unc_frac = 0.2, then the fetch length uncertainty will be 0.2 * fetch. Note that
    this is the fetch length uncertainty, NOT the uncertainty in Q due to fetch.
    """

    BASE   = 'w10'      # fixed wind‑speed variable name
    UNITS  = 'm_per_s'  # fixed unit specifier

    windspeed_primary_key       = f'{BASE}_{primary_windspeed_source}_{UNITS}'
    windspeed_secondary_key     = f'{BASE}_{secondary_windspeed_source}_{UNITS}'
    windspeed_std_primary_key   = f'{BASE}_std_{primary_windspeed_source}_{UNITS}'
    windspeed_std_secondary_key = f'{BASE}_std_{secondary_windspeed_source}_{UNITS}'
    windspeed_key               = f'{BASE}_{primary_windspeed_source}_{secondary_windspeed_source}_{UNITS}'
    windspeed_std_key           = f'{BASE}_std_{primary_windspeed_source}_{secondary_windspeed_source}_{UNITS}'

    df = pd.read_csv(simple_ime_csv_filename)
    dfw = pd.read_csv(wind_csv_filename)

    # Insert a datetime column derived from plume_id
    df['datetime_UTC'] = pd.to_datetime([x.split('_')[0][3:].replace('t', ' ') for x in df['plume_id']])

    # Replace NaNs in the C_Q column with 0
    nans = np.isnan(df[c_q_key])
    df.loc[nans, c_q_key] = 0

    dfw[windspeed_key]     = dfw[windspeed_primary_key].fillna(dfw[windspeed_secondary_key])
    dfw[windspeed_std_key] = dfw[windspeed_std_primary_key].fillna(dfw[windspeed_std_secondary_key])

    dfw['windspeed_source'] = np.where(dfw[windspeed_primary_key].notna(), primary_windspeed_source, secondary_windspeed_source)

    df = df.set_index('plume_id').join(dfw.set_index('plume_id'), how='left', rsuffix='_w')

    df = df[np.isfinite(df[windspeed_key])]

    Q, sigma_Q, sigma_C, sigma_w, sigma_f = \
        compute_Q_and_uncertainty_utils.compute_Q_and_uncertainty(
            df[c_q_key],
            df[sum_c2_key],
            df[windspeed_key],
            df[windspeed_std_key],
            df[fetch_key],
            fetch_unc_frac)

    df['Q_kg_hr']      = Q
    df['Q_unc_kg_hr']  = sigma_Q
    df['sigma_C']      = sigma_C
    df['sigma_f']      = sigma_f
    df['sigma_w']      = sigma_w

    if os.path.exists(output_csv_filename):
        raise ValueError(f'The specified output file {output_csv_filename} already exists.')
    df.to_csv(output_csv_filename)

if __name__ == '__main__':
    compute_Q_and_Qunc_simple_IME()