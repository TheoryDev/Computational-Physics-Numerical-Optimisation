# -*- coding: utf-8 -*-
"""
The script verifies the parabolic minmiser works on cosh(x)
and that the gradient method works on f(x,y)=()**2+()**2+5
"""

import optimiser
import numpy as np

f=optimiser.coshTest # cosh(x) to be minimised by parabolic minmiser
optimiser.parabolicTest(f=f,tol=0.000001)
# It returns a result very close to zero verifying the parabolic minmiser works.
g=optimiser.Test2D # f(x,y) to be minmised by gradient method
optimiser.dualOptimiserGradientTest(f=g)