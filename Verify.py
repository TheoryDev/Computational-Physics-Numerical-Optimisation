# -*- coding: utf-8 -*-
"""
The script calcualtes the mean error, it fits the fitfunction with no background to a normalised histogram
of the 10,000 lifetime measurements. It reproduces figure 1
"""

import dataRead
import fit
import optimiser
import numpy as np

w=fit.mData()
meanSig=w.meanSigma()
w.plotFitFunctionTest(sigma=meanSig) #reproduces figure 1
w.normalisationTest(sigma=meanSig)
