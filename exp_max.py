"""
Functions for performing the expectation maximisation algorithm over the
observed spike-pattern rates and the natural parameters. These functions
use 'call by reference' in that, rather than returning a result, they update
the data referred to by their parameters.

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

import probability
import transforms



CONVERGED = 1+1e-6



def compute_A(sigma_t0, sigma_t1, F):
    """
    TODO explain what A is

    :param numpy.ndarray sigma_t0:
        Covariance of the theta distribution for timestep t.
    :param numpy.ndarray sigma_t1:
        Covariance of the theta distribution for timestep t+1.
    :param numpy.ndarray F:
        Autoregressive parameter of state transisions, of dimensions (D, D).

    :returns:
        TODO explain what A is
    """
    a = numpy.dot(sigma_t0, F.T)
    A = numpy.dot(a, numpy.linalg.inv(sigma_t1))

    return A


def e_step(emd):
    """
    Computes the expectation (a multivariate Gaussian distribution) of the
    natural parameters of observed spike patterns, given the state-transition
    hyperparameters. Firstly performs a `forward' iteration, in which the
    expectation at time t is determined from the observed patterns at time t and
    the expectation at time t-1. Secondly performs a `backward' iteration, in
    which these sequential expectation estimates are smoothed over time.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    """
    # Compute the 'forward' filter density
    e_step_filter(emd)
    # Compute the 'backward' smooth density
    e_step_smooth(emd)


def e_step_filter(emd):
    """
    Computes the one-step-prediction density and the filter density in the
    expectation step.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    """
    # Iterate forwards over each timestep, computing filter density
    emd.theta_f[0,:], emd.sigma_f[0,:] = emd.max_posterior(emd, 0)
    for i in xrange(1, emd.T):
        # Compute one-step prediction density
        emd.theta_o[i,:] = numpy.dot(emd.F, emd.theta_f[i-1,:])
        tmp = numpy.dot(emd.F, emd.sigma_f[i-1,:,:])
        emd.sigma_o[i,:,:] = numpy.dot(tmp, emd.F.T) + emd.Q
        # Get MAP estimate of filter density
        emd.theta_f[i,:], emd.sigma_f[i,:] = emd.max_posterior(emd, i)


def e_step_smooth(emd):
    """
    Computes smooth density in the expectation step.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    """
    # Initialise the smoothed theta and sigma values
    emd.theta_s[-1,:] = emd.theta_f[-1,:]
    emd.sigma_s[-1,:,:] = emd.sigma_f[-1,:,:]
    # Iterate backwards over each timestep, computing smooth density
    for i in reversed(xrange(emd.T - 1)):
        # Compute the A matrix
        A = compute_A(emd.sigma_f[i,:,:], emd.sigma_o[i+1,:,:], emd.F)
        # Compute the backward-smoothed means
        tmp = numpy.dot(A, emd.theta_s[i+1,:] - emd.theta_o[i+1,:])
        emd.theta_s[i,:] = emd.theta_f[i,:] + tmp
        # Compute the backward-smoothed covariances
        tmp = numpy.dot(A, emd.sigma_s[i+1,:,:] - emd.sigma_o[i+1,:,:])
        tmp = numpy.dot(tmp, A.T)
        emd.sigma_s[i,:,:] = emd.sigma_f[i,:,:] + tmp


def m_step(emd):
    """
    Computes the optimised hyperparameters of the natural parameters of the
    expectation distributions over time. `Q' is the covariance matrix of the
    transition probability distribution. `F' is the autoregressive parameter of
    the state transitions, but it is kept constant in this implementation.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    """
    # Update the initial mean of the one-step-prediction density
    emd.theta_o[0,:] = emd.theta_s[0,:]
    # Compute the state-transition hyperparameter
    m_step_Q(emd)


def m_step_F(emd):
    """
    Computes the optimised autogregressive hyperparameter `F' of the natural
    parameters of the expectation distributions over time. See equation 39 of
    the source paper for details.

    NB: This function is not called in this implementation because the
    autoregressive parameter is kept constant.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    """
    # Set up temporary-results arrays
    a = numpy.zeros((emd.D, emd.D))
    b = numpy.zeros((emd.D, emd.D))
    # Sum partial results over each timestep
    for i in xrange(1, emd.T):
        A = compute_A(emd.sigma_s[i-1,:,:], emd.sigma_s[i,:,:], emd.F)
        a += numpy.dot(A, emd.sigma_s[i,:,:]) +\
             numpy.outer(emd.theta_s[i,:], emd.theta_s[i-1,:])
        b += emd.sigma_s[i-1,:,:] +\
             numpy.outer(emd.theta_s[i-1,:], emd.theta_s[i-1,:])
    # Dot the results
    emd.F = numpy.dot(a, numpy.linalg.inv(b))


def m_step_Q(emd):
    """
    Computes the optimised state-transition covariance hyperparameters `Q' of
    the natural parameters of the expectation distributions over time.

    :param container.EMData emd:
        All data pertaining to the EM algorithm.
    """
    lmbda = 0
    for i in xrange(1, emd.T):
        A = compute_A(emd.sigma_f[i-1,:,:], emd.sigma_o[i,:,:], emd.F)
        lag_one_covariance = numpy.dot(A, emd.sigma_s[i,:])
        lmbda += numpy.trace(emd.sigma_s[i,:,:]) +\
                 numpy.dot(emd.theta_s[i,:], emd.theta_s[i,:]) -\
                 2 * numpy.trace(lag_one_covariance) -\
                 2 * numpy.dot(emd.theta_s[i-1,:], emd.theta_s[i,:]) +\
                 numpy.trace(emd.sigma_s[i-1,:,:]) +\
                 numpy.dot(emd.theta_s[i-1,:], emd.theta_s[i-1,:])
    emd.Q = lmbda / emd.D / (emd.T - 1) * numpy.identity(emd.D)
