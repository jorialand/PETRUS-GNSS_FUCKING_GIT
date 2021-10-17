
from collections import OrderedDict
from math import sqrt
from scipy.special import erfinv

def updateMin(CurrentMin, Value):
    return min(CurrentMin, Value)

def updateMax(CurrentMax, Value):
    return max(CurrentMax, Value)

def updateHist(Hist, Value, Resolution):
    Bin = float(int(Value/Resolution)) * Resolution
    if Bin in Hist:
        Hist[Bin] = Hist[Bin] + 1

    else:
        Hist[Bin] = 1

def computeCdfFromHistogram(Histogram, NSamples):
    # Sort histogram
    SortedHist = OrderedDict({})
    for key in sorted(Histogram.keys()):
        SortedHist[key] = Histogram[key]
    # Define internal variables
    Cdf = OrderedDict({})
    Sigmas = OrderedDict({})
    CumulatedSamples = 0
    for Bin, Samples in SortedHist.items():
        CumulatedSamples = CumulatedSamples + Samples
        # Compute Cdf
        Cdf[Bin] = float(CumulatedSamples) / NSamples
        # Compute Sigmas
        Sigmas[Bin] = float(Bin/(sqrt(2)*erfinv(Cdf[Bin])))

    return Cdf, Sigmas

def computePercentile(Cdf, Percentile):
    for Bin, Freq in Cdf.items():
        if (Freq * 100) > Percentile:
            return Bin

def computeOverbound(Sigmas, ThresholdBin):
    SigmaOver = 0.0
    for Bin, Sigma in Sigmas.items():
        if Bin >= ThresholdBin:
            SigmaOver = updateMax(SigmaOver, Sigma)

    return SigmaOver
