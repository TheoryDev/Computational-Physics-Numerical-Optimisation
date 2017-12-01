
"""
The fit Module contains the class mData. The class mData reads the data from the file lifetime.txt 
by importing the function to do so from dataRead.py.
The class computes the fit function in the cases of no background and where background is included.
The class produces a plot of the fit function is superimpose on a histogram of the lifetime measurements.
The fit function is integrated over the range of its domain in order to verify it is normalised.
The values of the NLL is calculated for particular values of the parameters in the cases of no background and where background is included
The NLL is plotted as a function of the paramaters in the cases of no background and where background is included.

"""
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as spi
import scipy.special
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

class mData:   
    """ The class mData allows
    """
    def __init__(self):        
        """ The class object is initialised with attributes: An array of the 10,000 lifetime measurements 
            and their associated errors, constants B and C and N are used in
            calculations manytimes.
        """
        self.C= 1./np.sqrt(2) # constants C and B are used to calculates 
        self.B=1./np.sqrt(2*np.pi)
        self.data = np.loadtxt("lifetime.txt")
        self.N=len(self.data)    
        
        
    def meanSigma(self):
        """calculates the mean errror of the 10,000 lifetime measurements"""
        mean = np.mean(self.data[:,1]) 
        print ('Mean Error= ', mean)
        return mean    
        
    def plotFitFunctionTest(self,start=-5,end=5,step=0.03,tau=0.4,sigma=0.3): 
        """ start,end,step determined the lower limit,upper limit and step size in the numpy array of the lifetime 
            measurements that are plotted. tau and sigma are the parameters of the fit function.
            The function superimposes a plot of the fitfunction on a normalised histogram of 
            the 10,000 lifetime measurements. 
        """
        t=np.arange(start,end,step) # list of lifetimes to be plotted by the fitfunction.         
        fm=self.fitFunction(t,sigma,tau)    
        plt.figure()
        plt.plot(t,fm,label='Theoretical Model') # plot of fitfunction over range t
        plt.hist(self.data[:,0],bins='sqrt',normed=1,label='Histogram of lifetimes') # Histogram of 10,000 lifetime measurements   
        plt.xlabel('t/ps')
        plt.ylabel('fm(t)')
        plt.title('Fit of theoretical model to histogram of lifetime')
        plt.legend()
        plt.show()
        
    def normalisationTest(self,tau=0.4,sigma=0.2): 
        """ tau is the mean lifetime, sigam is the error 
            Checks the integral of the pdf over all time is equal to unity
        """
        value= spi.quad(lambda t: self.fitFunction(t,sigma,tau),-50*sigma,np.inf) 
        print ('Integral of fit function over all time= ', value[0])
        # Integrating from -inf to inf causes problems due to numerical rounding issues.  
        return value[0] # To remedy this the lower limit is -50*sigma
        
    def fitFunction(self,t,sigma,tau): 
        """t is the lifetime, sigma is the error and tau is the mean lifetime.
           The function calcualates the fitfunction value for a 
           given lifetime with parameters the mean error and the mean lifetime    
        """
        B=sigma/tau       
        fmi=(0.5/tau)*scipy.exp(0.5*B*B-(t/tau))*scipy.special.erfc(self.C*(B-(t/sigma)))# fit function value
        return fmi
        
    def fitDataFunction(self,graph=False,tau=0.4): 
        """ tau is the mean life time, if graph==True a graph is plotted
            Calculates the fit function for each of the 10,000 lifetime measurements. """                          
        self.fm=self.fitFunction(self.data[:,0],self.data[:,1],tau) # calculates fit function        
        if graph==True:
            plt.plot(self.data[:,0],self.fm,'yo')
            plt.xlabel('t/ps')
            plt.ylabel('fm(t)')            
        return self.fm    
    

    def fitDataFunctionOpt(self,x,graph=False,N=10000): 
        """ x is an array/list of tau values. The fit function is calculated
            from the first to Nth lifetime measurement. The method returns
            an array of the NLLs for each tau value in x.
        """
        points=[] # will store NLL for each tau
        for i in x:     
            self.fm=self.fitFunction(self.data[:N,0],self.data[:N,1],tau=i) # calculates array of fitfunction values            
            Log =-np.log(self.fm)    
            points.append(np.sum(Log)) #NLL value is appended.       
        return np.array(points)# array of NLL values
       
    def NLL(self,tau0=0.1,tauN=5,num=1000): 
        """tau0,tauN and num are used to create a linspace array where the start is tau0
           the end is tauN and the num is num. The method produces a plot of the NLL as a function of 
           tau.
        """
        self.nLogs=[]         
        taus = np.linspace(tau0,tauN,num)  
        for i in taus:             
            nLog = -np.log(self.fitDataFunction(tau=i))    
            self.nLogs.append(np.sum(nLog)) 
        plt.figure()
        plt.plot(taus,self.nLogs)
        plt.xlabel('tau/ps')
        plt.ylabel('NLL')
        plt.title('NLL as a function of tau')
        
    def fitFunctionBackground(self,t,sigma,tau,a=0.5):    
        """ t is the lifetime, sigma is the error and tau is the mean lifetime.
            a is the fraction of the signal in the sample.
            The method returns the fitfunction value in the case where background is included
        """
        A=t/sigma        
        back= self.B*(np.exp(-0.5*A*A)/sigma)# intermediate value from gaussian background
        fmib= a*self.fitFunction(t,sigma,tau)+(1.-a)*back # fitfunction with background      
        return fmib

        
    def NLLContourPlotBackground(self,tau0=0.39,tauN=0.42,num=100,a0=0.95,aN=1.0,num2=100,contourStep=0.05,N=20,Figure4=False,Figure5= False): 
        """ tau0,tauN and num are used to create a linspace array of taus where the start is tau0
            the end is tauN and the num is num.  and a0,aN and num2 are used to create a linspaces of a's where
            a0 is start, aN is end , num1 is num. The method produces a contour plot of the NLL as a function of 
            a and tau.
            If Figure4 == True, Figure 4 from the report is plotted.
            If Figure5 == True, Figure 5 from the report is plotted.
            """
        if Figure4 == True and Figure5 == True:
            raise Exception('Cannot plot both the small and large graphs at the same time')
        if Figure4 == True:    
            tau0,tauN,num,a0,aN,num2,N,contourStep=0.2,2,100,0.1,1.0,100,20,0.2
        if Figure5 == True:
            tau0,tauN,num,a0,aN,num2,N,contourStep=0.390,0.420,100,0.95,1.00,100,20,0.005
      
        nLogs=[]         
        taus = np.linspace(tau0,tauN,num)
        aS= np.linspace(a0,aN,num2)               
        for i in taus: 
            tmp=[]         
            for j in aS:                
                nLog = -np.log(self.fitFunctionBackground(self.data[:,0],self.data[:,1],tau=i,a=j))              
                tmp.append(np.sum(nLog))           
            nLogs.append(tmp)      
        X,Y = np.meshgrid(taus,aS) # meshgrid of taus and aS
      
        Z=np.array(nLogs)
        Z=np.transpose(Z)   # need transpose for it to work with meshgrid
        
        plt.figure()
        
        CS=plt.contour(X,Y,Z,N) # contour plot
        plt.clabel(CS,fontsize=12, color='k')
        plt.xlabel('tau/ps')
        plt.ylabel('a')
        plt.title('Contour plot of NLL(tau,a)')
        plt.ylim([a0,aN])
        plt.xlim([tau0,tauN])
        plt.xticks(np.arange(tau0,tauN,contourStep))         
        plt.show()
   
    def NLLBackground(self,tau,a):
        """Calculates the NLL value for a specific a and tau.
           tau is the mean lifetime and a is the fraction of the sample in the background
        """
        nll= -np.log(self.fitFunctionBackground(self.data[:,0],self.data[:,1],tau,a))
        return np.sum(nll)
                                            