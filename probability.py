"""
Functions for computing probabilities and distributions of spike patterns and
hyperparameters.

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
import pdb

import transforms



def log_likelihood(y_t, theta_f_t, R):
    """
    Computes the likelihood of observed pattern rates given the natural
    parameters, all for a single timestep.

    :param numpy.ndarray y_t:
        Frequency of observed patterns for one timestep.
    :param numpy.ndarray theta_f_t:
        Natural parameters of observed patterns for one timestep.
    :param int R:
        Number of trials over which patterns were observed.

    :returns:
        Log likelhood of the observed patterns given the natural parameters,
        as a float.
    """
    psi = transforms.compute_psi(theta_f_t)
    log_p = R * (numpy.dot(y_t, theta_f_t) - psi)

    return log_p


def log_marginal(emd, period=None):
    """
    Computes the log marginal probability of the observed spike-pattern rates
    by marginalising over the natural-parameter distributions. See equation 45
    of the source paper for details.

    This is just a wrapper function for `log_marginal_raw`. It unpacks data
    from the EMD container pbject and calls that function.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    :param period tuple:
        Timestep range over which to compute probability.

    :returns:
        Log marginal probability of the synchrony estimate as a float.
    """
    # Unwrap the parameters and call the raw function
    log_p = log_marginal_raw(emd.theta_f, emd.theta_o, emd.sigma_f, emd.sigma_o,
        emd.y, emd.R, period)

    return log_p


def log_marginal_raw(theta_f, theta_o, sigma_f, sigma_o, y, R, period=None):
    """
    Computes the log marginal probability of the observed spike-pattern rates
    by marginalising over the natural-parameter distributions. See equation 45
    of the source paper for details.

    From within SSASC, this function should be accessed by calling
    `log_marginal` with the EMD container as a parameter. This raw function is
    designed to be called from outside SSASC, when a complete EMD container
    might not be available.

    See the container.py for a full description of the parameter properties.

    :param period tuple:
        Timestep range over which to compute probability.

    :returns:
        Log marginal probability of the synchrony estimate as a float.
    """
    if period == None: period = (0, theta_f.shape[0])
    # Initialise
    log_p = 0
    # Iterate over each timestep and compute...
    a, b = 0, 0
    for i in xrange(period[0], period[1]):
        a += log_likelihood(y[i,:], theta_f[i,:], R)
        theta_d = theta_f[i,:] - theta_o[i,:]
        sig_inv = numpy.linalg.inv(sigma_o[i,:,:])
        b -= numpy.dot(theta_d, numpy.dot(sig_inv, theta_d))
        b += numpy.log(numpy.linalg.det(sigma_f[i,:,:])) -\
             numpy.log(numpy.linalg.det(sigma_o[i,:,:]))
    log_p = a + b / 2

    return log_p


def log_multivariate_normal(x, mu, sigma):
    """
    Computes the probability of `x' from a multivariate normal distribution
    of mean `mu' and covariance `sigma'. This function is taken from SciPy.

    :param float x:
        Point at which to evalue the multivariate normal PDF.
    :param numpy.ndarray mu:
        Mean of the multivariate normal, of size D.
    :param numpy.ndarray sigma:
        Covariance of the multivariate normal, of dimensions (D, D).

    :returns:
        Probability of x as a float.
    """
    # Find the number of dimensions of the multivariate normal
    dim = mu.size
    assert mu.size == sigma.shape[0] == sigma.shape[1],\
        'The covariance matrix must have the same number of dimensions as '+\
        'the mean vector.'
    # Compute the normalising term
    norm = numpy.log(((2 * numpy.pi) ** dim * numpy.linalg.det(sigma)) ** .5)
    # Compute the distribution term
    si, xd = numpy.linalg.inv(sigma), x - mu
    dist = -.5 * numpy.dot(numpy.dot(xd.transpose(), si), xd)
    # Finish up!
    log_p = norm + dist

    return log_p
