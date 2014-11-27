"""
Functions for generating synthetic spike data.

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



def generate_spikes(p, R, seed=None):
    """
    Draws spike patterns for each of `R' trial runs from the probability mass
    specified in `p'. `p' must have T rows, one independent probability mass for
    each timestep, and 2^C columns, where C is the number of cells (maximum of
    8) involved in the spike pattern.

    :param numpy.ndarray p:
        Probability mass of spike patterns for each timestep.
    :param int R:
        Number of spike patterns to generate for each timestep.

    :returns:
        Binary matrix with dimensions (time, runs, cells), in which a `1' in
        location (t, r, c) denotes a spike at time t in run r by cell c, as a
        numpy.ndarray
    """
    # Set metadata
    T, N = p.shape[0], numpy.int(numpy.log2(p.shape[1]))
    # Initialise random seed
    numpy.random.seed(seed)
    # Set up spike patterns
    fx = transforms.enumerate_patterns(N)
    # Set up the output array (time, trials, cells)
    T, C = p.shape[0], numpy.log2(p.shape[1])
    spikes = numpy.zeros((T, R, N))
    # Iterate over each probability
    for i in xrange(T):
        # Draw random values from the probability distribution
        idx = random_weighted(p[i], R)
        # Extract spike patterns for each trial
        spikes[i,:,:] = fx[idx,:]

    return spikes


def random_weighted(p, R):
    """
    Draws `R' integers from the probability mass over the integers `p'.

    :param numpy.ndarray p:
        Probability mass.
    :param int R:
        Sample size to draw from `p'.

    :returns:
        `R' random numbers drawn from distribution `p', as a numpy.ndarray.
    """
    # Take a cumulative sum of the probability mass
    cs = numpy.cumsum(p)
    # Draw uniform random numbers for each timestep or trial
    rnd = numpy.random.random(R)
    # For each random value, find the index of the first weight above it
    idx = numpy.zeros(R, dtype=numpy.int)
    for i in xrange(R):
        idx[i] = numpy.sum(cs < rnd[i])

    return idx
