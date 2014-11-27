"""
Minimal working example of the SSASC program, examining second-order interaction
between two cells. Note that modules are imported at the top of each section
that uses them (in contrast to the usual convention of importing all modules at
the top of the file) so as to be completely explicit about the external
requirements.

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

# Set time (milliseconds), number of trials, and number of cells
T, R, N = 1000, 100, 2
# Set the interaction order
O = 2


# ----- SPIKE SYNTHESIS -----
# Global module
import numpy
# Local modules
import synthesis
import transforms

# Create theta parameters for a single timestep
theta = numpy.array([-3, -3, 3]) # i.e. [theta_1, theta_2, theta_12]
# Repeat the parameters for every timestep
theta = numpy.tile(theta, T).reshape(T, theta.size)
# Initialise the transforms library in preparation for computing P
transforms.initialise(N, O)
# Compute P for each timestep
p = numpy.zeros((T, 2**N))
for i in xrange(T):
    p[i,:] = transforms.compute_p(theta[i,:])
# Generate spikes!
spikes = synthesis.generate_spikes(p, R, seed=1)


# ----- ALGORITHM EXECUTION -----
# Global module
import numpy
# Local module
import __init__ # From outside this folder, this would be 'import ssasc'

# Run the algorithm!
emd = __init__.run(spikes, O, lmbda=.005)


# ----- PLOTTING -----
# Global module
import pylab

# Set up an output figure
fig, ax = pylab.subplots(2, 1, sharex=1)
# Plot theta traces
ax[0].plot(emd.theta_s[:,0], c='b')
ax[0].plot(emd.theta_s[:,1], c='r')
ax[1].plot(emd.theta_s[:,2], c='g')
# Set axis limits
ax[0].set_ylim(-4.5, -1.5)
ax[1].set_ylim(1.5, 4.5)
# Set labels
ax[0].set_title('Second order interaction between two cells')
ax[0].set_ylabel('First-order theta')
ax[1].set_xlabel('Time (ms)')
ax[1].set_ylabel('Second-order theta')
# Show figure!
pylab.show()
