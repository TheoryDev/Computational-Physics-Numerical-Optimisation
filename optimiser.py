# -*- coding: utf-8 -*-
"""
This module is used to generate the results. It contains methods 
to implement a parabolic one dimensional minimum search. It contains a method
to implement a two dimensional minimum search.
It calculates errors in 1d using the last parabolic estimate or changing the NLL by 0.5 units.
It calculates errors in 2d by changing the NLL by 0.5 units
It also test the optimisers on suitable functions
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import fit
from scipy.optimize import curve_fit
import scipy as sp

p=math.pow(10,-12)
class optimise(fit.mData):
    """The optimise class inherits from the fit.mData class.
       It then inherits from the fit It contains methods 
       to implement a parabolic one dimensional minimum search. It contains a method
       to implement a two dimensional minimum search.
       It calculates errors in 1d using the last parabolic estimate or changing the NLL by 0.5 units.
       It calculates errors in 2d by changing the NLL by 0.5 units

    """
    def __init__(self):        
        """The class object is instantiated with attributes: An array containing the 10,000 lifetime measurements 
            and their respective errors. The optimised atribute is used by methods to prevent unecessary calculatuons.
            C and B that are frequently used constants. """
        self.C= 1./np.sqrt(2) # make global C
        self.data = np.loadtxt("lifetime.txt")
        self.B=1./np.sqrt(2*np.pi)
        self.N=1  
        self.optimised= False
    
        
    def optim(self,x0=[0.1,0.3,0.8],tol=0.00001,N=10000,Print=True): #parabolic optimiser
        """ This method takes arguments: x0-the initial contions, tol-tolerance level, N-the number of measurements in the NLL fit """
        x=np.array(x0)    
        xDif=tol*10.0 
        tmp = 0.
        while abs(xDif) > tol:  # repeats the loop until the difference of two successive iterations of x3 is within the tolerence level.
            x=np.sort(x)
            y=self.fitDataFunctionOpt(x,N=N)             
            tmp1= (x[2]*x[2]-x[1]*x[1])*y[0]+(x[0]*x[0]-x[2]*x[2])*y[1]+(x[1]*x[1]-x[0]*x[0])*y[2]
            tmp2= (x[2]-x[1])*y[0]+(x[0]-x[2])*y[1]+(x[1]-x[0])*y[2]
            x3 = tmp1/(2.0*tmp2) # calculates x3
            xDif=tmp-x3 # difference between current and previous x3's
            tmp = x3        
            y3 = self.fitDataFunctionOpt(np.array([x3]),N=N)
            self.tmpX=x
            x,y = np.append(x,x3),np.append(y,y3[0])
            inMax= np.argmax(y)
            x,y= np.delete(x,inMax),np.delete(y,inMax) # removes maximum of 4 points
            self.x=x
        self.minTau=x[np.argmin(y)] 
        if Print == True:
            print ('Tau at the minimum= ' , self.minTau)
        return self.minTau # returns the minimum
        
    def accuracyHalfNLL(self,x0=[0.1,0.3,0.8],h=0.00001,tol=0.000000001):
        """x0 is the initial conditions, h is the step size, tol is the tolerence level
           The method calcualtes the values of tau at which the NLL changes by 0.5 units away from the minimum.
           It then returns the errors.
        """
        minTau = self.optim(x0,tol,Print=False)
        nll= self.fitDataFunctionOpt(np.array([minTau]))
      
        tauPlus=minTau+h*minTau 
        tauMinus=minTau-h*minTau
        
        while 1==1:        
            nllPlus=self.fitDataFunctionOpt(np.array([tauPlus]))
            nllDifP = nllPlus[0] - nll[0]       
            if 0.5< abs(nllDifP) :               
                break
            tauPlus=tauPlus+h*minTau #creates next tau value to try    
       
        while 1==1:        
            nllMinus=self.fitDataFunctionOpt(np.array([tauMinus]))
            nllDifM = nllMinus[0] - nll[0]      
            if 0.5< abs(nllDifM) :     
                break
            tauMinus=tauMinus-h*minTau #creates next tau value to try
      
        sigmaP,sigmaM = minTau-tauMinus, tauPlus-minTau
        print ('Positive Error in Tau at Minimum=', sigmaP)
        print ('Negative Error in Tau at Minimum=', sigmaM)
        
        return sigmaP,sigmaM
    
    def accuracyCurvature(self,x0=[0.1,0.3,0.8],tol=0.00000001,N=10000,Print=True):
        """ x0 is the initial conditions, tol is the tolerence level, N is the number of measurements
            included in the NLL fit. The method returns the 
        """
       # print N
        minTau = self.optim(x0,tol,N=N,Print=Print)  # gets the last parabolic estimate and it becomes a object attribute
        y=self.fitDataFunctionOpt(self.tmpX,N=N)   
        x=self.tmpX
        d=(x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1])
        curv=(2./d)*(y[0]*(x[2]-x[1])+y[1]*(x[0]-x[2])+y[2]*(x[1]-x[0]))#calculates second derivative
        self.error=1./math.sqrt(abs(curv)) # estimates error from curvature
        if Print==True:
            print ('Curvature Error in Tau at Minimum= ' , self.error)
        return self.error

    def rootN(self,x,a):
        """ x is a input number, a is a constant of proportionality. Do not confused with fraction of signal in sample  """
        return a/np.sqrt(x)
    
    def nPoints(self,x0=[0.1,0.3,0.8],tol=0.000001):
        """x0 is the initial coordinates,tol is the tolerence level. 
           The method produces a plot of the error as a function of N and fits a curve to it.
           It returns the value of N needed for an error of 10^-15 seconds by extrapolation.
        """
        plt.figure()
        Ns=np.arange(5000,10000,100)                 
        errors=[]      
        for i in Ns: # calculates the error for element in Ns      
            errors.append(self.accuracyCurvature(x0,tol,N=i, Print=False)) 
        # i is the number of measurments included into the NLL fit
        self.popt, self.pcov = curve_fit(self.rootN,Ns,errors)# optimises a/sqrt(N) plot to fit data  
        nMax= math.pow((self.popt[0]/(0.001)),2) # N needed for error of 10^-15 seconds
        N2=np.arange(5000,nMax+500,100)
        rootNs=1./np.sqrt(N2)
        errorPlot =plt.plot(Ns,errors,'bx',label='errors') # plot of errorrs vs N
        rootPlot= plt.plot(N2,self.popt[0]*rootNs,'g-',label=' N^-0.5 curve fit') # plot of optimised curve       
        plt.xlabel('N')
        plt.ylabel('error/ps')
        plt.title('error vs N')
        plt.legend()
        print ('Number of points needed for error of 10^-15s=' , nMax)
        plt.show()
        return nMax
        
    def dualOptimiserGradient(self,aG=0.985,tauG=0.410,h=0.0001,alpha=0.00001,graph=True):
        """ aG is the initial guess of a,tauG is the initial guess of tau, h is the step size, 
            alpha is a constant in the gradient method, if graph == True, a graph of the minimum
            search is plotted.
        """
        self.optimised=True # prevents other methods calling this method if it has beeen done already.
       
        a=aG
        tau=tauG    
        tmp=1.+h
        counter=0
        nllDif = 10.
        TauP,Ap,nllP =[],[],[]
        while nllDif >= 0.:
            counter = counter+1
            nll= self.NLLBackground(tau,a) # current value of NLL
            dndtau=(self.NLLBackground(tau*tmp,a)-nll)/h # partial derivative of NLL with respect to tau
            dnda=(self.NLLBackground(tau,a*tmp)-nll)/h # partial derivative of NLL with respect to a
           
     
            tmp1,tmp2=tau,a 
            tau=tau-alpha*dndtau # gradient method for tau
            a=a-alpha*dnda # gradient method for a           
            tmp3=self.NLLBackground(tau,a) # new value of NLL
        
            nllDif=nll-tmp3
            TauP.append(tmp1)
            Ap.append(tmp2)
            nllP.append(nll)
        Ap.append(a)
        nllP.append(tmp3)
        TauP.append(tau)
        if graph==True:
            # prints the minimum search        
            plt.figure('search for minimum NLL')
            plt.plot(TauP,Ap,'-.')
            plt.xlabel('tau')
            plt.ylabel('a')    
        print ('Tau at Minimum= ',tmp1)
        print ('Fraction of background at the Minimum= ', 1- tmp2)
        self.tauMin=tmp1
        self.aMin=tmp2
        return tmp1, 1-tmp2    
   
        
    def ErrorHalfNLLBackground(self,h=0.000003):  
        """     """        
        if self.optimised == False:            
            self.dualOptimiserGradient(graph=False,h=h)  
        tauMin= self.tauMin
        aMin= self.aMin
        tauPlus=tauMin+h*tauMin
        
        nll = self.NLLBackground(tau=tauMin,a=aMin)
        while 1==1:        
            nllPlus=self.NLLBackground(tau=tauPlus,a=aMin)
            nllDifP = nllPlus - nll   
            #print nllDifP
            if 0.5 < nllDifP :            
                break
            tauPlus=tauPlus+h*tauMin     
        tauMinus=tauMin-h*tauMin
        while 1==1:        
            nllMinus=self.NLLBackground(tau=tauMinus,a=aMin)
            nllDifM = nllMinus - nll      
            if 0.5< nllDifM  :    
                break
            tauMinus=tauMinus-h*tauMin
      
        sigmaPTau,sigmaMTau = tauMin-tauMinus, tauPlus-tauMin     
        
        print ('Positive Error in Tau at Minimum=', sigmaPTau)
        print ('Negative Error in Tau at Minimum=', sigmaMTau)
        nllDifP,nllDifM = 0.,0.
        aPlus=aMin+h*aMin
        while 1==1:        
            nllPlus=self.NLLBackground(tau=tauMin,a=aPlus)
            nllDifP = nllPlus - nll             
            if 0.5 < nllDifP :            
                break
            aPlus=aPlus+h*aMin     
        aMinus=aMin-h*aMin
        while 1==1:        
            nllMinus=self.NLLBackground(tau=tauMin,a=aMinus)
            nllDifM = nllMinus - nll      
            if 0.5< nllDifM  :    
                break
            aMinus=aMinus-h*aMin
        sigmaPa,sigmaMa= aPlus-aMin, aMin-aMinus
        # It follows on from the propagation of errors that the error in 1-a is equal to the error in a.
        print ('Positive Error in (1-a) at Minimum=',sigmaPa)
        print ('Negative Error in (1-a) Minimum=', sigmaMa)          
       
    
def dualOptimiserGradientTest(f,x0=1.5,y0=1.5, h=0.00000001,alpha=0.001,graph=True):
    """f is the function to be minimised, x0 and y0 are the initial guesses of x and y,
       alpha is the constant in the gradient method, if graph== True a plot of the minimum search is created.
       The method tests the gradient method algorithm on a suitable function to verify it.       
    """
    x=x0
    y=y0    
    tmp=1.+h   
    
    uDif=0.001
    xS,yS =[],[]
    while uDif>0.: # test has converged when the next iteration increases the function's value
        
        u=f(x,y)        
        dudx=(f(x*tmp,y)-u)/h # partial derivative of function with respect to x
        dudy=(f(x,y*tmp)-u)/h # partial derivative of function with respect to y
 
        tmpx,tmpy=x,y        
        x=x-alpha*dudx # gradient method on x
        y=y-alpha*dudy # gradient method on y
        tmp3=f(x,y)
        
        uDif=u-tmp3 # difference in the value of the function at sucessive iterations
        xS.append(x)
        yS.append(y)
      
    if graph==True:          
        plt.figure('search for minimum u(x,y)')
        plt.plot(xS,yS,'-.')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    print ('x at minimum= ', tmpx)
    print ('y at minimum= ', tmpy)
    print ('function at minimum= ', u)
    return tmpx,tmpy   
    
def parabolicTest(f,x0=[0.1,0.3,0.8],tol=0.00001):
    """ This method takes arguments: f the function to be minimised, x0-the initial contions, tol-tolerance level.
        It returns the minimum of the function
    """
    x=np.array(x0)    
    xDif=tol*10.0 
    tmp = 0.
    while abs(xDif) > tol:  # repeats the loop until the difference of two successive iterations of x3 is within the tolerence level.
        x=np.sort(x)
        y=f(x)             
        tmp1= (x[2]*x[2]-x[1]*x[1])*y[0]+(x[0]*x[0]-x[2]*x[2])*y[1]+(x[1]*x[1]-x[0]*x[0])*y[2]
        tmp2= (x[2]-x[1])*y[0]+(x[0]-x[2])*y[1]+(x[1]-x[0])*y[2]
        x3 = tmp1/(2.0*tmp2) # calculates x3
        xDif=tmp-x3 # difference between current and previous x3's
        tmp = x3        
        y3 = f(x3)
        tmpX=x
        x,y = np.append(x,x3),np.append(y,y3)
        inMax= np.argmax(y)
        x,y= np.delete(x,inMax),np.delete(y,inMax) # removes maximum of 4 points
        
    minX= x[np.argmin(y)] 
    print ('x at the minimum= ', minX)
    print ('function at the minimum= ' , np.argmin(y))
    return minX # returns the minimum

def Test2D(x,y):
    "Takes x and y values are returns the value of the function"
    return (x-1.)*(x-1.)+(y-1.)*(y-1.)+5.        
                    
def coshTest(x):
    """" Takes x values and returns cosh x """
    return np.cosh(x)
            
        
    
        