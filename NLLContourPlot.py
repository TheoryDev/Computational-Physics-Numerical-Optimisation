# -*- coding: utf-8 -*-
"""
The script returns figures 4 and Figure 5
"""

import dataRead
import fit
import optimiser
import numpy as np


w=fit.mData()
w.NLLContourPlotBackground(Figure4=True, Figure5=False)# If Figure 4==True, Figure 4 is reproduced, 
#If figure 5==True Figure 5 is reproduced. 
