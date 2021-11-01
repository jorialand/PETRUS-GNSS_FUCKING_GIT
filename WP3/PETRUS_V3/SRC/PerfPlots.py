#!/usr/bin/env python

import sys, os
from pandas import read_csv
from pandas import DataFrame

sys.path.append(os.getcwd() + '/' + \
                os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
from COMMON.Plots import generatePlot
from InputOutput import HistIdx, PerfIdx
import numpy as np
from scipy import stats
from collections import OrderedDict

# Performances plots configuration flags
ConfPerf = OrderedDict({})
ConfPerf["PLOT_VPE_HISTOGRAM"] = 1#
ConfPerf["PLOT_AVAILABILITY"] = 1#
ConfPerf["PLOT_CONT_RISK"] = 1#
ConfPerf["PLOT_HPE95"] = 1#
ConfPerf["PLOT_VPE95"] = 1#
ConfPerf["PLOT_MAX_HSI"] = 1#
ConfPerf["PLOT_MAX_VSI"] = 1#
ConfPerf["PLOT_MIN_HPL"] = 1#
ConfPerf["PLOT_MIN_VPL"] = 1#
ConfPerf["PLOT_MAX_HPL"] = 1#
ConfPerf["PLOT_MAX_VPL"] = 1#
ConfPerf["PLOT_MIN_SATNUM"] = 1#
ConfPerf["PLOT_MAX_SATNUM"] = 1#
ConfPerf["PLOT_MAX_HDOP"] = 1#
ConfPerf["PLOT_MAX_VDOP"] = 1#

def initPlot(PerfFilesList, PlotConf, Title, Label, Service):
    # Compute information from PerfFilesList
    PerfFileName = os.path.basename(PerfFilesList[0])
    PerfFileNameSplit = PerfFileName.split('_')
    DatepDat = PerfFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    # Dump information into PlotConf
    PlotConf["Title"] = "%s %s on Year %s DoY %s" % (Service, Title, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PERF/Figures/%s/' % Label + '%s_%s_Y%sD%s.png' % (Label, Service, Year, Doy)


def initHist(HistFile, PlotConf, Title, Label):
    # Compute information from HistFile
    HistFileName = os.path.basename(HistFile)
    HistFileNameSplit = HistFileName.split('_')
    Rcvr = HistFileNameSplit[2]
    DatepDat = HistFileNameSplit[3]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    # Dump information into PlotConf
    PlotConf["xLabel"] = "LPV200 VPE Histogram"

    PlotConf["Title"] = "%s %s on Year %s DoY %s" % (Rcvr, Title, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PERF/Figures/%s/' % Label + '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)


# Plot availability map
def plotAvailability(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Availability Percentage", "AVAIL", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Avail = PerfData[PerfIdx["AVAIL"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Availability Percentage [%]"
    PlotConf["ColorBarMin"] = min(Avail)
    PlotConf["ColorBarMax"] = 100.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["nData"] = {}
    PlotConf["xData"][0] = Lon
    PlotConf["yData"][0] = Lat
    PlotConf["zData"][0] = Avail
    PlotConf["nData"][0] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot continuity risk map
def plotContRisk(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Continuity Risk", "CONT_RISK", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["Decimal"] = True

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Cont = PerfData[PerfIdx["CONTRISK"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Continuity Risk"
    PlotConf["ColorBarMin"] = min(Cont)
    PlotConf["ColorBarMax"] = max(Cont)
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Cont
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot HPE 95% map
def plotHPE95(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "HPE 95%", "HPE_95%", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Hpe95 = PerfData[PerfIdx["HPE95"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "HPE 95% [m]"
    PlotConf["ColorBarMin"] = min(Hpe95)
    PlotConf["ColorBarMax"] = max(Hpe95)
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Hpe95
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot VPE 95% map
def plotVPE95(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "VPE 95%", "VPE_95%", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Vpe95 = PerfData[PerfIdx["VPE95"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "VPE 95% [m]"
    PlotConf["ColorBarMin"] = min(Vpe95)
    PlotConf["ColorBarMax"] = max(Vpe95)
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Vpe95
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot extrapolated VPE map
def plotExtVPE(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Extrapolated VPE", "EXT_VPE", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    ExtVpe = PerfData[PerfIdx["EXTVPE"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Extrapolated VPE [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 10.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = ExtVpe
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum HSI map
def plotMaxHSI(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum HSI", "MAX_HSI", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Hsi = PerfData[PerfIdx["HSIMAX"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum HSI [m]"
    PlotConf["ColorBarMin"] = min(Hsi)
    PlotConf["ColorBarMax"] = max(Hsi)
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Hsi
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum VSI map
def plotMaxVSI(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum VSI", "MAX_VSI", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Vsi = PerfData[PerfIdx["VSIMAX"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum VSI [m]"
    PlotConf["ColorBarMin"] = min(Vsi)
    PlotConf["ColorBarMax"] = max(Vsi)
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Vsi
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot minimum HPL map
def plotMinHPL(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Minimum HPL", "MIN_HPL", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Hpl = PerfData[PerfIdx["HPLMIN"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Minimum HPL [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 20.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Hpl
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot minimum VPL map
def plotMinVPL(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Minimum VPL", "MIN_VPL", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Vpl = PerfData[PerfIdx["VPLMIN"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Minimum VPL [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 20.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Vpl
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum HPL map
def plotMaxHPL(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum HPL", "MAX_HPL", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Hpl = PerfData[PerfIdx["HPLMAX"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum HPL [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 200.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Hpl
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum VPL map
def plotMaxVPL(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum VPL", "MAX_VPL", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Vpl = PerfData[PerfIdx["VPLMAX"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum VPL [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 200.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Vpl
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot minimum number of satellites map
def plotMinSats(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Minimum Number of Satellites", "MIN_SATSNUM", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["Integer"] = True

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    MinSats = PerfData[PerfIdx["NSVMIN"]][FilterCond].to_numpy(dtype=np.float64)
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Minimum Number of Satellites"
    PlotConf["ColorBarMin"] = min(MinSats)
    PlotConf["ColorBarMax"] = max(MinSats)
    PlotConf["ColorBarTicks"] = range(int(min(MinSats)), int(max(MinSats)) + 1)

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = MinSats
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum number of satellites map
def plotMaxSats(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum Number of Satellites", "MAX_SATSNUM", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["Integer"] = True

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    MaxSats = PerfData[PerfIdx["NSVMAX"]][FilterCond].to_numpy(dtype=np.float64)
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum Number of Satellites"
    PlotConf["ColorBarMin"] = min(MaxSats)
    PlotConf["ColorBarMax"] = max(MaxSats)
    PlotConf["ColorBarTicks"] = range(int(min(MaxSats)), int(max(MaxSats)) + 1)

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = MaxSats
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum HDOP map
def plotMaxHDOP(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum HDOP", "MAX_HDOP", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Hdop = PerfData[PerfIdx["HDOPMAX"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum HDOP [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 20.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Hdop
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot maximum VDOP map
def plotMaxVDOP(Service, PerfFilesList, PerfData):
    # Graph settings definition
    PlotConf = {}
    initPlot(PerfFilesList, PlotConf, "Maximum VDOP", "MAX_VDOP", Service)

    PlotConf["Type"] = "Map"
    PlotConf["FigSize"] = (12.6, 10.4)

    PlotConf["LonMin"] = -35
    PlotConf["LonMax"] = 40
    PlotConf["LatMin"] = 10
    PlotConf["LatMax"] = 85
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = 'o'

    # Prepare data to be plotted
    FilterCond = PerfData[PerfIdx["SERVICE"]] == Service
    Lon = PerfData[PerfIdx["LON"]][FilterCond].to_numpy()
    Lat = PerfData[PerfIdx["LAT"]][FilterCond].to_numpy()
    Vdop = PerfData[PerfIdx["VDOPMAX"]][FilterCond].to_numpy()
    Rcvr = PerfData[PerfIdx["RCVR"]][FilterCond].to_numpy()

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Maximum VDOP [m]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 20.0
    # PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = Lon
    PlotConf["yData"] = Lat
    PlotConf["zData"] = Vdop
    PlotConf["nData"] = Rcvr

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# PerfPlots main functions
# ----------------------------------------------------------

def generateHistPlot(LPV200ExtVpe, HistFile):
    # LPV200 VPE HISTOGRAM
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_VPE_HISTOGRAM"] == 1):
        # Read the cols we need from HistFile file
        HistData = read_csv(HistFile, delim_whitespace=True, skiprows=1, header=None, \
                            usecols=[HistIdx["BINMIN"], HistIdx["BINMAX"], HistIdx["BINFREQ"]])

        print('Plot LPV200 VPE Histogram ...')

        # Graph settings definition
        PlotConf = {}
        initHist(HistFile, PlotConf, "LPV200 VPE Histogram", "LPV200_VPE_HISTOGRAM")

        PlotConf["Type"] = "Hist"
        PlotConf["FigSize"] = (14.4, 8.4)

        PlotConf["yLabel"] = "Relative Frequency PDF"

        PlotConf["Grid"] = True
        PlotConf["BarWidth"] = GnssConstants.HIST_RES

        # Prepare data to be plotted
        BinMin = HistData[HistIdx["BINMIN"]].to_numpy()
        BinMax = HistData[HistIdx["BINMAX"]].to_numpy()
        HistData = HistData[HistIdx["BINFREQ"]].to_numpy()

        # Compute area under the histogram
        AreaHist = 0.0
        for element in HistData:
            AreaHist = AreaHist + GnssConstants.HIST_RES * element

        # Generate overbounding curve
        xOver = np.arange(-max(BinMax), max(BinMax), GnssConstants.HIST_RES)
        yOver = stats.norm.pdf(xOver, 0, LPV200ExtVpe / 5.33)

        # Normalize overbounding curve
        Norm = AreaHist
        xOverbound = []
        yOverbound = []
        for i in range(len(xOver)):
            if xOver[i] >= 0:
                xOverbound.append(xOver[i])
                yOverbound.append(2 * Norm * yOver[i])

        PlotConf["Legend"] = ["Gaussian Overbounding", "Min: " + str(min(BinMin)) + "\n" + "Max: " + str(max(BinMax))]

        # Plotting histogram
        PlotConf["xData"] = BinMin
        PlotConf["yData"] = HistData

        # Plotting overbounding curve
        PlotConf["x2Data"] = xOverbound
        PlotConf["y2Data"] = yOverbound

        # Call generatePlot from Plots library
        generatePlot(PlotConf)


def generatePerfPlots(Service, PerfFilesList):
    # AVAILABILITY MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_AVAILABILITY"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["AVAIL"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["AVAIL"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Availability Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotAvailability(Service, PerfFilesList, PerfData)

    # CONTINUITY RISK MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_CONT_RISK"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["CONTRISK"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["CONTRISK"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Continuity Risk Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotContRisk(Service, PerfFilesList, PerfData)

    # HPE 95% MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_HPE95"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["HPE95"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["HPE95"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot HPE 95% Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotHPE95(Service, PerfFilesList, PerfData)

    # VPE 95% MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_VPE95"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["VPE95"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["VPE95"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot VPE 95% Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotVPE95(Service, PerfFilesList, PerfData)

    # MAXIMUM HSI MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_HSI"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["HSIMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["HSIMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum HSI Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxHSI(Service, PerfFilesList, PerfData)

    # MAXIMUM VSI MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_VSI"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["VSIMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["VSIMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum VSI Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxVSI(Service, PerfFilesList, PerfData)

    # MINIMUM HPL MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MIN_HPL"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["HPLMIN"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["HPLMIN"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Minimum HPL Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMinHPL(Service, PerfFilesList, PerfData)

    # MINIMUM VPL MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MIN_VPL"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["VPLMIN"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["VPLMIN"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Minimum VPL Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMinVPL(Service, PerfFilesList, PerfData)

    # MAXIMUM HPL MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_HPL"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["HPLMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["HPLMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum HPL Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxHPL(Service, PerfFilesList, PerfData)

    # MAXIMUM VPL MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_VPL"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["VPLMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["VPLMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum VPL Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxVPL(Service, PerfFilesList, PerfData)

    # MINIMUM NUMBER OF SATELLITES
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MIN_SATNUM"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["NSVMIN"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["NSVMIN"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Minimum Number of Satellites Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMinSats(Service, PerfFilesList, PerfData)

    # MAXIMUM NUMBER OF SATELLITES
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_SATNUM"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["NSVMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["NSVMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum Number of Satellites Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxSats(Service, PerfFilesList, PerfData)

    # MAXIMUM HDOP MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_HDOP"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["HDOPMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["HDOPMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum HDOP Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxHDOP(Service, PerfFilesList, PerfData)

    # MAXIMUM VDOP MAP
    # ----------------------------------------------------------
    if (ConfPerf["PLOT_MAX_VDOP"] == 1):
        # Initialize PerfData Dataframe
        PerfData = DataFrame(
            columns=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"], PerfIdx["VDOPMAX"]])
        # Read the cols we need from all PerfFile files
        for PerfFile in PerfFilesList:
            NewData = read_csv(PerfFile, delim_whitespace=True, skiprows=1, header=None, \
                               usecols=[PerfIdx["RCVR"], PerfIdx["LON"], PerfIdx["LAT"], PerfIdx["SERVICE"],
                                        PerfIdx["VDOPMAX"]])
            # Append information to PerfData
            PerfData = PerfData.append(NewData, ignore_index=True)

        print('Plot Maximum VDOP Map in ' + Service + '...')

        # Configure plot and call plot generation function
        plotMaxVDOP(Service, PerfFilesList, PerfData)
