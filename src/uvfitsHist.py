# This routine takes a MWA uvfits file produced by cotter and evaluates
# the coverage of the fractional floating point space. This can be used
# to estimate the potential for lossless compression, not taking into
# account digital noise introduced by the algorithm.
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import argparse

def ind_array(refarray, comparray):
    """
    Return an array of the indices of every element of comparray
    in refarray.
    """
    rsorted = np.argsort(refarray)
    cpos = np.searchsorted(refarray[rsorted], comparray)
    ind = rsorted[cpos]
    return ind 

def fextract( f ):
    """
    Function extracts mantissa, exponent and sign from floating point 
    number f.
    """
    bits = np.asarray( f, dtype=np.float64 ).view( np.int64 )
    if not bits & 0x7fffffffffffffff: # f == +/-0
        return 0, 0
    sign = np.sign(bits)
    exponent = ( (bits>>52) & 0x7ff ) - 1075
    mantissa = 0x10000000000000 | ( bits & 0xfffffffffffff )
    # from here on f == sign * mantissa * 2**exponent
    for shift in 32, 16, 8, 4, 2, 1:
        if not mantissa & ((1<<shift)-1):
            mantissa >>= shift
            exponent += shift
    return sign * mantissa, exponent

def uvfitsSplit(ff, N):
    """
    Function to split a uvfits file into N parts, where N <= GCOUNT.
    ff is a pyfits instance.
    
    Returns N or N+1 pyfits instances. If float(GCOUNT)/N - GCOUNT/N <= 0.5
    the last part will contain more, else the last part will contain less
    random groups.
    """
    nparts = ff[0].header['GCOUNT']/N
    if (ff[0].header['GCOUNT']/float(N) - nparts)> 0.5:
        N += 1
    narr = np.ones(N,'i32') * nparts
    narr[-1] = ff[0].header['GCOUNT'] - (N-1)*nparts
    #*** There must be a bug here!!!***
    for i in narr:
        pyfits.open()
    return narr

def main_loop(data, cnt=1000):
    """
    Main worker function
    """
    tlen = 0
    for i in range(cnt):
        real=data[i][5][0][0][:][:,0,0]
        imag=data[i][5][0][0][:][:,0,1]
        ampl = np.sqrt(real**2+imag**2)
        val = np.abs(ampl)# we may want do this separately for real and imag
        tlen += len(val)
        frac = val - np.int32(val)
        ufrac = np.unique(frac)
        if i == 0:
            utfrac = ufrac #init utfrac array
            ind = ind_array(utfrac, frac)
        else:
            nfracs = np.lib.arraysetops.setdiff1d(ufrac,utfrac) # new fractions
            utfrac = np.append(utfrac,nfracs) # concatenate fracs
            ind_n = ind_array(utfrac, frac) # get all indices
            ind = np.append(ind, ind_n) # add indexes
    
    return utfrac, ind, np.bincount(ind), tlen

def main_loop_orig(data, cnt=1000):
    for i in range(cnt):
        real=data[i][5][0][0][:][:,0,0]
        imag=data[i][5][0][0][:][:,0,1]
        ampl = np.sqrt(real**2+imag**2)
        val = np.abs(ampl)
        frac = val - np.int32(val)
        frac.sort()
        if i ==0:
            tfrac = frac
        else:
            tfrac = np.append(frac,tfrac)


    utfrac = np.unique(tfrac)
    dig = np.digitize(tfrac,utfrac,right=True) # fraction indexes
    dig = np.bincount(dig)  # histogram of fractions appearance
    dig = np.bincount(dig)  # histogram of reappearance of fractions
    return utfrac, dig, len(tfrac)

def coveragePlot(bins):
    """
    Plots the histogram of the values found in the uvfits file.
    """
    plt.bar(range(len(bins)),bins)
    return

if __name__ == "__main__":
    #These should be input paramaters:
    fname = '/Users/awicenec/tmp/1125780432.uvfits'

    parser = argparse.ArgumentParser(description="Histogram of uvfits data fractions")
    parser.add_argument("-f", "--fileName", default=fname, help="The name of the uvfits file", type=str)
    parser.add_argument("-c", "--count", help="Number of groups to read", default=100, type=int)
    a = parser.parse_args()

 
    ff=pyfits.open(a.fileName)
    data = ff[0].data
    utfrac, ind, hist1, tlen = main_loop(data, cnt=a.count)
    hist2 = np.bincount(hist1)
    print "Ratio all/unique fractions: %4.2f" % (float(tlen)/len(utfrac))
    coveragePlot(hist2[:np.where(hist2[1:]==0)[0][0]])
    plt.show()
