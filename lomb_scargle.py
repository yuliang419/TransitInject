import numpy as np
from astroML.time_series import search_frequencies, lomb_scargle, MultiTermFit


# Compute the best frequencies
def compute_best_frequencies(data, n_eval=10000, n_retry=5, generalized=True):
    """
    Find best peaks with LS periodogram
    :param data: [time, flux, uncertainty]
    :param n_eval: The number of points to evaluate in the range on each iteration.
    :param n_retry: Number of top points to search on each iteration.
    :param generalized:
    :return:
    [omega, power]
    """
    t, y, dy = data
    kwargs = dict(generalized=generalized)
    omega, power = search_frequencies(t, y, dy, n_eval=n_eval,
                                      n_retry=n_retry,
                                      LS_kwargs=kwargs)
    results = [omega, power]

    return results


def fit_curve(data, results):
    """
    Reconstruct time series with the first few Fourier components
    :param data: [time, flux, uncertainty]
    :param results: [omega, power]
    :return:
    """

    t, y, dy = data
    omega, power = results
    omega_best = omega[np.argmax(power)]
    print " - omega_0 = %.10g" % omega_best

    # do a fit to the first 4 Fourier components
    mtf = MultiTermFit(omega_best, 4)
    mtf.fit(t, y, dy)
    phase_fit, y_fit = mtf.predict(1000)