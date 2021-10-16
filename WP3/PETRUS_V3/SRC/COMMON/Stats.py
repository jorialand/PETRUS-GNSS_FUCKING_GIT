from collections import OrderedDict
from numpy import sqrt
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
    Cdf = OrderedDict({})
    CumulatedSamples = 0
    for Bin, Samples in SortedHist.items():
        CumulatedSamples = CumulatedSamples + Samples
        Cdf[Bin] = float(CumulatedSamples) / NSamples

    return Cdf

def computePercentile(Cdf, Percentile):
    for Bin, Freq in Cdf.items():
        if (Freq * 100) > Percentile:
            return Bin

    return Cdf
