import numpy as np
import celerite
from celerite import terms
from scipy.optimize import minimize
from gatspy.periodic import LombScargleFast
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel


def get_period(time, flux):
    model = LombScargleFast(fit_offset=True).fit(time, flux)
    _, _ = model.periodogram_auto(oversampling=100)
    model.optimizer.period_range = (0.1, 30)
    print 'Best-fit period:', model.best_period
    return model.best_period


class CustomTerm(terms.Term):
    parameter_names = ("log_a", "log_b", "log_c", "log_P")

    def get_real_coefficients(self):
        b = np.exp(self.log_b)
        return (
            np.exp(self.log_a) * (1.0 + b) / (2.0 + b),
            np.exp(self.log_c),
        )

    def get_complex_coefficients(self):
        b = np.exp(self.log_b)
        return (
            np.exp(self.log_a) / (2.0 + b),
            0.0,
            np.exp(self.log_c),
            2*np.pi*np.exp(-self.log_P),
        )


def neg_log_like(params, y, gp, package='celerite'):
    if package == 'celerite':
        gp.set_parameter_vector(params)
        return -gp.log_likelihood(y)
    elif package == 'george':
        gp.kernel[:] = params
        ll = gp.lnlikelihood(y, quiet=True)
        # The scipy optimizer doesn't play well with infinities.
        return -ll if np.isfinite(ll) else 1e25
    else:
        raise ValueError('Invalid GP package name')


def celerite_fit(time, flux):
    # estimate period of periodic component
    period = get_period(time, flux)

    # now do GP
    bounds = dict(log_a=(None, None), log_b=(None, None), log_c=(None, -np.log(period)),
                  log_P=(np.log(period/2), None))
    sigma = np.std(flux[:50])
    sorted = np.sort(flux)
    amp = sorted[-3] - sorted[3]
    kernel = CustomTerm(log_a=np.log(amp), log_b=np.log(0.1), log_c=np.log(5), log_P=np.log(period),
                        bounds=bounds)
    gp = celerite.GP(kernel)
    gp.compute(time, sigma)
    print 'sigma=', sigma
    print("Initial log likelihood: {0}".format(gp.log_likelihood(flux)))

    initial_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()
    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(flux, gp))
    gp.set_parameter_vector(r.x)
    print r

    pred_mean, var = gp.predict(flux, time, return_var=True)
    sigma = np.std(flux-pred_mean)
    good = np.where(abs(flux-pred_mean) < 3*sigma)

    gp.compute(time[good], sigma)
    initial_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()
    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(flux[good], gp))
    gp.set_parameter_vector(r.x)
    params = gp.get_parameter_dict()
    print params

    print("Final log likelihood: {0}".format(gp.log_likelihood(flux[good])))
    pred_mean, var = gp.predict(flux[good], time, return_var=True)
    return pred_mean, var


def george_fit(time, flux):
    period = get_period(time, flux)

    sigma = np.std(flux[:50])
    sorted = np.sort(flux)
    amp = sorted[-3] - sorted[3]
    a, l, G, P = amp/2, 5, 0.1, period
    kernel = a * ExpSquaredKernel(l) * ExpSine2Kernel(G, P)
    gp = george.GP(kernel, mean=np.mean(flux))
    gp.compute(time, sigma)
    print 'Initial ln likelihood:', gp.lnlikelihood(flux)

    bounds = ((None, None), (np.log(period), None), (None, None), (np.log(period/2), None))
    initial_params = gp.kernel.vector
    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(flux, gp, 'george'))
    print r
    gp.kernel[:] = r.x
    pred_mean, var = gp.predict(flux, time)

    sigma = np.std(flux - pred_mean)
    good = np.where(abs(flux - pred_mean) < 3 * sigma)
    gp.compute(time[good], sigma)
    initial_params = gp.kernel.vector
    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(flux[good], gp, 'george'))
    gp.kernel[:] = r.x
    pred_mean, var = gp.predict(flux[good], time)
    print 'Final ln likelihood:', gp.lnlikelihood(flux[good])

    return pred_mean, var