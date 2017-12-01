
"""
Module used to read the lifetime measurements from the file.
"""
import matplotlib.pyplot as plt
import numpy as np

def loadData(graphs=False):
    """ The loadData method reads the 10,000 lifetime measurements and their associated errrors from the file lifetime.txt
        It returns the data. If graphs == True, then a histogram of the lifetime measurements is generated.
    """
    data=np.loadtxt("lifetime.txt")
  
    if graphs == True:
        plt.figure()
        plt.hist(data[:,0],bins='sqrt') # creates histogram       
        plt.title('Histogram of lifetimes')
        plt.xlabel('t/ps')
        plt.ylabel('frequency of bin')
    return data