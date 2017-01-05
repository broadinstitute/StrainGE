#import pandas as pd
#import matplotlib.pyplot as plt

def makeCorrPlot(name, xmin, fullsims, hashsims, oldhashsims=None):
    print "new r =", fullsims.corr(hashsims)
    if type(oldhashsims) != type(None):
        plt.scatter(fullsims, oldhashsims, color='red', label='Old')
        print "new r =", fullsims.corr(oldhashsims)
    plt.scatter(fullsims, hashsims, color='blue', label='New')
    plt.plot([0,1], color='black')
    plt.xlim(xmin, 1)
    plt.ylim(xmin, 1)
    plt.xlabel("Full Similarity Score")
    plt.ylabel("Fingerprint Similarity Score")
    plt.legend(loc=4)
    plt.savefig(name)
    plt.show()

def makeErrorPlot(name, xmin, fullsims, hashsims, oldhashsims=None):
    if type(oldhashsims) != type(None):
        plt.scatter(fullsims, oldhashsims - fullsims, color='red', label='Old')
    plt.scatter(fullsims, hashsims - fullsims, color='blue', label='New')
    plt.xlim(xmin, 1)
    plt.plot([0,0], color='black')
    plt.xlabel("Full Similarity Score")
    plt.ylabel("Fingerprint Similarity Score Error")
    plt.legend(loc=4)
    plt.savefig(name)
    plt.show()

def makePlots(fullsims, hashsims, oldhashsims = None, cutoff=0.9):
    makeCorrPlot("corr.png", .2, fullsims, hashsims, oldhashsims)
    makeErrorPlot("error.png", .2, fullsims, hashsims, oldhashsims)
    f90 = fullsims > cutoff
    makeCorrPlot("corr90.png", cutoff, fullsims[f90], hashsims[f90], oldhashsims[f90])
    makeErrorPlot("error90.png", cutoff, fullsims[f90], hashsims[f90], oldhashsims[f90])





