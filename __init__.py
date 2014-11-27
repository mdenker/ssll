"""
Master file of the State-Space Analysis of Spike Correlations.

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

import container
import exp_max
import probability
import max_posterior
import synthesis
import transforms



def run(spikes, order, window=1, map_function='nr', lmbda=0.005, max_iter=10):
    """
    Master-function of the State-Space Analysis of Spike Correlation package.
    Uses the expectation-maximisation algorithm to find the probability
    distributions of natural parameters of spike-train interactions over time.
    Calls slave functions to perform the expectation and maximisation steps
    repeatedly until the data likelihood reaches an asymptotic value.

    Note that the execution of some slave functions to this master function are
    of exponential complexity with respect to the `order' parameter.

    :param numpy.ndarray spikes:
        Binary matrix with dimensions (time, runs, cells), in which a `1' in
        location (t, r, c) denotes a spike at time t in run r by cell c.
    :param int order:
        Order of spike-train interactions to estimate, for example, 2 =
        pairwise, 3 = triplet-wise...
    :param int window:
        Bin-width for counting spikes, in milliseconds.
    :param string map_function:
        Name of the function to use for maximum a-posterior estimation of the
        natural parameters at each timestep. Refer to max_posterior.py.
    :param float lmdbda:
        Coefficient on the identity matrix of the initial state-transition
        covariance matrix.
    :param int max_iter:
        Maximum number of iterations for which to run the EM algorithm.

    :returns:
        Results encapsulated in a container.EMData object, containing the
        smoothed posterior probability distributions of the natural parameters
        of the spike-train interactions at each timestep, conditional upon the
        given spikes.
    """
    # Ensure NaNs are caught
    numpy.seterr(invalid='raise', under='raise')
    # Initialise the EM-data container
    map_func = max_posterior.functions[map_function]
    emd = container.EMData(spikes, order, window, map_func, lmbda)
    # Initialise the coordinate-transform maps
    transforms.initialise(emd.N, emd.order)
    # Set up loop guards for the EM algorithm
    lmp = -numpy.inf
    lmc = probability.log_marginal(emd)
    # Iterate the EM algorithm until convergence or failure
    while (emd.iterations < max_iter) and (emd.convergence > exp_max.CONVERGED):
        # Perform EM
        exp_max.e_step(emd)
        exp_max.m_step(emd)
        # Update previous and current log marginal values
        lmp = lmc
        lmc = probability.log_marginal(emd)
        # Update EM algorithm metadata
        emd.iterations += 1
        emd.convergence = lmp / lmc

    return emd
