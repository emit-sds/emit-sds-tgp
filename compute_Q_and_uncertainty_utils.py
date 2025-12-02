import numpy as np

def compute_Q_and_uncertainty(C_Q, sumC2, windspeed_m_s, windspeed_std_m_s, fetch_m, fetch_unc_fraction):
    """
    Compute the emission rate Q and its uncertainties. Uses standard error propagation of:
    Q = windspeed / fetch * sum_over_pixels(k * concentration_length)

    Parameters
    ----------
    C_Q : ndarray of floats
        Emission rate divided by windspeed [kg/hr / (m/s)]
    sumC2 : ndarray of floats
        Sum of squared concentration uncertainties across all pixels.
    windspeed_m_s : ndarray of floats
        Mean wind speed [m/s]
    windspeed_std_m_s : ndarray of floats
        Standard deviation of wind speed measurements in [m/s]
    fetch_m : ndarray of floats
        Fetch length [m]
    fetch_unc_fraction : ndarray of floats
        Fractional uncertainty in fetch length (e.g., 0.1 for 10 % uncertainty) [unitless]

    Returns
    -------
    Q : ndarray of floats
        Computed emission rate in kg/hr
    sigma_Q : ndarray of floats
        Standard deviation of total emission rate uncertainty in kg/hr
    sigma_C : ndarray of floats [kg/hr]
        Standard deviation of emission rate uncertainty due to concentration.
    sigma_w : ndarray of floats [kg/hr]
        Standard deviation of emission rate uncertainty due to wind speed.
    sigma_f : ndarray of floats [kg/hr]
        Standard deviation of emission rate uncertainty due to fetch length.
    """
    # Compute emission rate Q as conversion factor times wind speed
    Q = C_Q * windspeed_m_s

    # Pre‑compute term (C_Q * fetch_m)^2 for later uncertainty calculations
    C2 = (C_Q * fetch_m)**2

    # Uncertainty in fetch length (absolute value)
    fetch_std = fetch_unc_fraction * fetch_m

    # Squared uncertainty in Q due to concentration length uncertainty in each pixel
    sigma2_C = windspeed_m_s**2 / fetch_m**2 * sumC2

    # Squared uncertainty in Q due to fetch length uncertainty
    sigma2_f = windspeed_m_s**2 / fetch_m**4 * C2 * fetch_std**2

    # Squared uncertainty in Q due to wind speed uncertainty
    sigma2_w = 1 / fetch_m**2 * C2 * (windspeed_std_m_s)**2

    # Total squared uncertainty in Q (sum of independent contributions)
    sigma2_Q = sigma2_C + sigma2_f + sigma2_w

    # Return emission rate and standard deviations (square roots of variance terms)
    return Q, np.sqrt(sigma2_Q), np.sqrt(sigma2_C), np.sqrt(sigma2_w), np.sqrt(sigma2_f)