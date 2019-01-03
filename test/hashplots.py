#  Copyright (c) 2016-2019, Broad Institute, Inc. All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither the name Broad Institute, Inc. nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

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





