# -*- coding: utf-8 -*-
"""
This script plots the NLL with no background it reproduces Figure 2
"""

import dataRead
import fit
import optimiser
import numpy as np

w=fit.mData()
w.NLL() # reproduces figure 2
