# -*- coding: utf-8 -*-
"""
The script produces the results in the case where background is included.
It calculates the values of tau and the fraction of the background in the 
signal at the minimum of the NLL. It then calculates their respective errors
by method of changing the NLL by 0.5 units.

"""

import dataRead
import fit
import optimiser
import numpy as np

a=optimiser.optimise()
a.dualOptimiserGradient()# finds minimum tau and a.
a.ErrorHalfNLLBackground()# finds errors in tau and a by chaning NLL by 0.5 units