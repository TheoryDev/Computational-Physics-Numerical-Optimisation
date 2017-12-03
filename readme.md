# In this project I estimated the mean life time of the a meson (subatomic particle) by minimising the negative log likelihood fit of 10,000 lifetime measurements.

The lifetime is first estimated with no background signal using a 1D parabolic minimiser.
Then the lifetime is estimated with a background signal using the minumum gradient descent method.

Modules used to model physics and statistics:
* dataRead.py
* fit.py
* optimiser.py

Data:
* lifetime.txt - contains 10,000 lifetime measurements and their associated errors.

Scripts - scripts are used to generate results and graphs:
* Histograms.py
* Verify.py
* NLLplotNoBackground.py
* ResultsNoBackground.py
* TestMinimisers.py
* NLLContourPlot.py
* ResultsWithBackground.py

## `dataRead.py:`

The loadData method reads the 10,000 lifetime measurements 
and their associated errrors from the file lifetime.txt. 
It returns the data in an array. 
If `graphs` == True, then a histogram of the lifetime measurements is generated.

## `fit.py:`

The fit Module contains the class mData. The class mData reads the data from the file lifetime.txt 
by importing the function to do so from dataRead.py.
The class computes the fit function in the cases of no background and where background is included.
The class produces a plot of the fit function is superimpose on a histogram of the lifetime measurements.
The fit function is integrated over the range of its domain in order to verify it is normalised.
The values of the NLL is calculated for particular values of the parameters 
in the cases of no background and where background is included
The NLL is plotted as a function of the paramaters 
in the cases of no background and where background is included.

## `optimiser.py:`

This module is used to generate the results. It contains methods 
to implement a parabolic one dimensional minimum search. It contains a method
to implement a two dimensional minimum search.
It calculates errors in 1d using the last parabolic estimate or changing the NLL by 0.5 units.
It calculates errors in 2d by changing the NLL by 0.5 units
It also can test the optimisers on suitable functions

## `Histograms.py:`

This script reads the 10,000 measurements data from the file lifetime.txt.
It produces histograms of the data

## `Verify.py:`

The script calcualtes the mean error, it fits the fitfunction with no background to a normalised histogram
of the 10,000 lifetime measurements, ensures the fitfunction is normalised and produces a plot.

<p align="center"> 
 <img src="/images/modelHisto.png" height= "500" width="500">
 </p>

## `NLLplotNoBackground.py:`

The produces a plot of the NLL as a function of tau.

<p align="center"> 
 <img src="/images/Figure1.png" height= "400" width="400">
 </p>

## `ResultsNoBackground.py:`

The script produces the results in the case where the
background is negleceted. It calculates the value of tau at
the minimum and the associated errors obtained by changing the 
NLL by 0.5 units and using the curvature of the last Lagrange polynomial.

 <p align="center"> 
 <img src="/images/NoBackground.png" height= "400" width="400">
 </p>

## `TestMinimisers.py:`

The script verifies the `parabolic minmiser` works on `cosh(x)`
and that the `gradient method` works on `f(x,y)=(x-1)**2+(y-1)**2+5`


## `NLLContourPlot.py:`

This script produces contour plots of the NLL as a 
funcition to tau and a. If the argument `Figure4` == True
Figure 4 is reproduced. If the argument `Figure5` == True
Figure 5 is reproduced. The argument corresponding to the
Figure you do not want must be equal to False or an 
exception will be raised.

 <p align="center"> 
 <img src="/images/contour.png" height= "400" width="400">
 </p>

## `ResultsWithBackground.py:`

The script produces the results in the case where background is included.
It calculates the values of tau and the fraction of the background in the 
signal at the minimum of the NLL. It then calculates their respective errors
by method of changing the NLL by 0.5 units.

 <p align="center"> 
 <img src="/images/search_for_minimum_NLL.png" height= "400" width="400">
 </p>
 
# `(note any figures refered to are in the context of the report for this project)`
