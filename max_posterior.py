"""
Functions for computing maximum a-posterior probability estimates of natural
parameters given the observed data.

---

State-Space Analysis of Spike Correlations (Shimazaki et al. PLoS Comp Bio 2012)
Copyright (C) 2014  Thomas Sharp (thomas.sharp@riken.jp)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy

import probability
import transforms



# Named function pointers to MAP estimators
# SEE BOTTOM OF FILE

# Parameters for gradient-ascent methods of MAP estimation
MAX_GA_ITERATIONS = 100
GA_CONVERGENCE = 1e-3



def newton_raphson(emd, t):
    """
    Computes the MAP estimate of the natural parameters at some timestep, given
    the observed spike patterns at that timestep and the one-step-prediction
    mean and covariance for the same timestep.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    :param int t:
        Timestep for which to compute the maximum posterior probability.

    :returns:
        Tuple containing the mean and covariance of the posterior probability
        density, each as a numpy.ndarray.
    """
    # Extract observed patterns and one-step predictions for time t
    y_t = emd.y[t,:]
    theta_o = emd.theta_o[t,:]
    sigma_o = emd.sigma_o[t,:,:]
    R = emd.R
    # Initialise theta_max to the smooth theta value of the previous iteration
    theta_max = emd.theta_s[t,:]
    # Use non-normalised posterior prob. as loop guard
    lpp = -numpy.inf
    lpc = probability.log_likelihood(y_t, theta_max, R) +\
          probability.log_multivariate_normal(theta_max, theta_o, sigma_o)
    iterations = 0
    # Iterate the gradient ascent algorithm until convergence or failure
    while lpc - lpp > GA_CONVERGENCE:
        # Compute the eta of the current theta values
        p = transforms.compute_p(theta_max)
        eta = transforms.compute_eta(p)
        # Compute the inverse of one-step covariance
        sigma_o_i = numpy.linalg.inv(sigma_o)
        # Compute the first derivative of the posterior prob. w.r.t. theta_max
        dllk = R * (y_t - eta)
        dlpr = -numpy.dot(sigma_o_i, theta_max - theta_o)
        dlpo = dllk + dlpr
        # Compute the second derivative of the posterior prob. w.r.t. theta_max
        ddlpo = -R * transforms.compute_fisher_info(p, eta) - sigma_o_i
        # Dot the results to climb the gradient, and accumulate the result
        ddlpo_i = numpy.linalg.inv(ddlpo)
        theta_max -= numpy.dot(ddlpo_i, dlpo)
        # Update previous and current posterior prob.
        lpp = lpc
        lpc = probability.log_likelihood(y_t, theta_max, R) +\
              probability.log_multivariate_normal(theta_max, theta_o, sigma_o)
        # Count iterations
        iterations += 1
        # Check for check for overrun
        if iterations == MAX_GA_ITERATIONS:
            raise Exception('The maximum-a-posterior gradient-ascent '+\
                'algorithm did not converge before reaching the maximum '+\
                'number iterations.')

    return theta_max, -ddlpo_i



# Named function pointers to MAP estimators
functions = {'nr': newton_raphson}
