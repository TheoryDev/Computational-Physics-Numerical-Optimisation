# -*- coding: utf-8 -*-
"""
The script produces the results in the case where the
background is negleceted. It calculates the value of tau at
the minimum and the associated errors obtained by changing the 
NLL by 0.5 units and using the curvature of the last Lagrange polynomial.
It also produces a plot of the error as a function of the number of measurements in the NLL fit N.
A curve is fitted to the plot and it is extrapolated to estimate the number of points needed
for an error of 10^-15 seconds that value is then returned.
It reproduces figure 2
"""

import dataRead
import fit
import optimiser
import numpy as np


a=optimiser.optimise()
a.optim() # finds minimim
a.accuracyCurvature()# finds error from curvature
a.accuracyHalfNLL() # finds errors by changing NLL by 0.5 units
a.nPoints() #reproduces figure 2 and calculates the Number of points need for an error of 10^-15 seconds
