import sys, os
from pandas import read_csv
from InputOutput import PosIdx

sys.path.append(os.getcwd() + '/' + \
                os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
from COMMON.Plots import generatePlot
import numpy as np
from scipy.stats import gaussian_kde
from collections import OrderedDict

# Activate/Deactivate plots
# Position plots configuration flags
ConfPos = OrderedDict({})
ConfPos["PLOT_DOP"] = 1#
ConfPos["PLOT_ERR_vs_LIM"] = 1#
ConfPos["PLOT_ERROR"] = 1#
ConfPos["PLOT_HPE_vs_HDOP"] = 1#
ConfPos["PLOT_SAF_INDEX"] = 1#
ConfPos["PLOT_HOR_STANDFORD"] = 1#
ConfPos["PLOT_VER_STANDFORD"] = 1#

def initPlot(PosFile, PlotConf, Title, Label, Sol):
    # Compute information from PosFile
    PosFileName = os.path.basename(PosFile)
    PosFileNameSplit = PosFileName.split('_')
    Rcvr = PosFileNameSplit[1]
    DatepDat = PosFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    # Dump information into PlotConf
    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s in %s from %s on Year %s" \
                        " DoY %s" % (Title, Sol, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/SPVT/Figures/%s/' % Label + \
                       '%s_%s_%s_Y%sD%s.png' % (Label, Sol, Rcvr, Year, Doy)


# Plot DOPS
def plotDops(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "DOPS", "DOPS_vs_TIME", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4, 7.6)
    PlotConf["SecondAxis"] = ["Num SV " + Sol[0]]

    PlotConf["yLabel"] = "DOP [m]"
    PlotConf["y2Label"] = "Number of satellites used in " + Sol[0] + " solution"

    PlotConf["y2Lim"] = [0, max(sorted(PosData[PosIdx["NVS-SOL"]])) + 1]

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = ["HDOP", "VDOP", "PDOP", "TDOP", "Num SV " + Sol[0]]

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    for Label in PlotConf["Legend"]:
        if Label == "Num SV " + Sol[0]:
            PlotConf["xData"][Label] = PosData[PosIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
            PlotConf["yData"][Label] = PosData[PosIdx["NVS-SOL"]][FilterCond]
        else:
            PlotConf["xData"][Label] = PosData[PosIdx["SOD"]] / GnssConstants.S_IN_H
            PlotConf["yData"][Label] = PosData[PosIdx[Label]]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Position Errors vs Position Limits
def plotErrorVsLimit(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "PE vs PL", "POS_ERROR_vs_LIMIT", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4, 7.6)

    PlotConf["yLabel"] = "Value [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = ["VPE", "HPE", "HPL", "VPL"]

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    for Label in PlotConf["Legend"]:
        PlotConf["xData"][Label] = PosData[PosIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = abs(PosData[PosIdx[Label]][FilterCond].to_numpy())

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Position Errors
def plotErrors(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "Position Errors", "POS_ERROR_vs_TIME", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4, 7.6)

    PlotConf["yLabel"] = "Position Errors [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = ["VPE", "HPE"]

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    for Label in PlotConf["Legend"]:
        PlotConf["xData"][Label] = PosData[PosIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = PosData[PosIdx[Label]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot HPE vs HDOP
def plotHorizontalPE(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "HPE vs DOP", "HPE_vs_HDOP", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4, 7.6)

    PlotConf["xLabel"] = "EPE [m]"
    PlotConf["yLabel"] = "NPE [m]"

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '+'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "HDOP"
    PlotConf["ColorBarMin"] = min(sorted(PosData[PosIdx["HDOP"]]))
    PlotConf["ColorBarMax"] = max(sorted(PosData[PosIdx["HDOP"]]))
    PlotConf["ColorBarTicks"] = None
    # PlotConf['ColorBarTicks'] = range(int(PlotConf['ColorBarMin']), int(PlotConf['ColorBarMin']) + 1)
    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    Label = 0
    PlotConf["xData"][Label] = PosData[PosIdx["EPE"]][FilterCond]
    PlotConf["yData"][Label] = PosData[PosIdx["NPE"]][FilterCond]
    PlotConf["zData"][Label] = PosData[PosIdx["HDOP"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Safety Index
def plotSafeIndex(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "Safety Index", "SAFE_INDEX_vs_TIME", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4, 7.6)

    PlotConf["yLabel"] = "Safety Index"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = True
    PlotConf["HLine"] = [(1.0, 0, 24)]
    PlotConf["Legend"] = ["VSI", "HSI"]

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    for Label in PlotConf["Legend"]:
        PlotConf["xData"][Label] = PosData[PosIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = abs(PosData[PosIdx[Label]][FilterCond].to_numpy())

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Horizontal Stanford Diagram
def plotHorStand(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "Horizontal Stanford Diagram", "HOR_STANDFORD_DIAGRAM", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4, 7.6)

    PlotConf["xLabel"] = "HPE [m]"
    PlotConf["yLabel"] = "HPL [m]"

    if Sol[0] == "PA":
        PlotConf["xLim"] = [0, 50]
        PlotConf["yLim"] = [0, 50]
        PlotConf["HLine"] = [(GnssConstants.LPV200_HAL, 0, 50), (GnssConstants.APVI_HAL, 0, 50)]
        PlotConf["VLine"] = [(GnssConstants.LPV200_HAL, 0, 50), (GnssConstants.APVI_HAL, 0, 50)]
        PlotConf["SLine"] = [np.array([0, 50]), np.array([0, 50])]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    # Processing data to be plotted
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    Hpe = PosData[PosIdx["HPE"]][FilterCond].to_numpy()
    Hpl = PosData[PosIdx["HPL"]][FilterCond].to_numpy()
    # Get point density data
    Den = np.vstack([Hpe, Hpl])
    ZData = gaussian_kde(Den)(Den)
    # Sort data by point density values
    idx = ZData.argsort()
    Hpe, Hpl, ZData = Hpe[idx], Hpl[idx], ZData[idx]

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Point Density (Number of Samples : " + str(len(Hpe)) + ")"
    PlotConf["ColorBarMin"] = min(ZData)
    PlotConf["ColorBarMax"] = max(ZData)
    PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = Hpe
    PlotConf["yData"][Label] = Hpl
    PlotConf["zData"][Label] = ZData

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Vertical Stanford Diagram
def plotVerStand(PosFile, PosData, Sol):
    # Graph settings definition
    PlotConf = {}
    initPlot(PosFile, PlotConf, "Vertical Stanford Diagram", "VER_STANDFORD_DIAGRAM", Sol[0])

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4, 7.6)

    PlotConf["xLabel"] = "VPE [m]"
    PlotConf["yLabel"] = "VPL [m]"

    if Sol[0] == "PA":
        PlotConf["xLim"] = [0, 60]
        PlotConf["yLim"] = [0, 60]
        PlotConf["HLine"] = [(GnssConstants.LPV200_VAL, 0, 60), (GnssConstants.APVI_VAL, 0, 60)]
        PlotConf["VLine"] = [(GnssConstants.LPV200_VAL, 0, 60), (GnssConstants.APVI_VAL, 0, 60)]
        PlotConf["SLine"] = [np.array([0, 60]), np.array([0, 60])]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    # Processing data to be plotted
    FilterCond = PosData[PosIdx["SOL"]] == Sol[1]
    Vpe = PosData[PosIdx["VPE"]][FilterCond].to_numpy()
    Vpl = PosData[PosIdx["VPL"]][FilterCond].to_numpy()
    # Get point density data
    Den = np.vstack([Vpe, Vpl])
    ZData = gaussian_kde(Den)(Den)
    # Sort data by point density values
    idx = ZData.argsort()
    Vpe, Vpl, ZData = Vpe[idx], Vpl[idx], ZData[idx]

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Point Density (Number of Samples : " + str(len(Vpe)) + ")"
    PlotConf["ColorBarMin"] = min(ZData)
    PlotConf["ColorBarMax"] = max(ZData)
    PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = abs(Vpe)
    PlotConf["yData"][Label] = Vpl
    PlotConf["zData"][Label] = ZData

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


def generatePosPlots(PosFile):
    Sol = ['PA', 1]

    # ----------------------------------------------------------
    # PLOTTING FUNCTIONS
    # ----------------------------------------------------------

    # DOPS vs TIME
    # ----------------------------------------------------------
    if (ConfPos["PLOT_DOP"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["SOD"], PosIdx["NVS-SOL"], PosIdx["HDOP"], PosIdx["VDOP"], PosIdx["PDOP"],
                                    PosIdx["TDOP"], PosIdx["SOL"]])

        print('Plot DOPS vs Time in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotDops(PosFile, PosData, Sol)

    # POSITION ERRORS VS POSITION LIMITS
    # ----------------------------------------------------------
    if (ConfPos["PLOT_ERR_vs_LIM"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["SOD"], PosIdx["HPE"], PosIdx["VPE"], PosIdx["HPL"], PosIdx["VPL"],
                                    PosIdx["SOL"]])

        print('Plot Position Errors vs Position Limits vs Time in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotErrorVsLimit(PosFile, PosData, Sol)

    # POSITION ERRORS vs TIME
    # ----------------------------------------------------------
    if (ConfPos["PLOT_ERROR"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["SOD"], PosIdx["HPE"], PosIdx["VPE"], PosIdx["SOL"]])

        print('Plot Position Errors vs Time in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotErrors(PosFile, PosData, Sol)

    # HORIZONTAL POSITION ERROR vs HDOP
    # ----------------------------------------------------------
    if (ConfPos["PLOT_HPE_vs_HDOP"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["EPE"], PosIdx["NPE"], PosIdx["HDOP"], PosIdx["SOL"]])

        print('Plot Horizontal Position Error vs HDOP in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotHorizontalPE(PosFile, PosData, Sol)

    # SAFETY INDEX vs TIME
    # ----------------------------------------------------------
    if (ConfPos["PLOT_SAF_INDEX"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["SOD"], PosIdx["HSI"], PosIdx["VSI"], PosIdx["SOL"]])

        print('Plot Safety Index vs Time in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotSafeIndex(PosFile, PosData, Sol)

    # HORIZONTAL Stanford DIAGRAM
    # ----------------------------------------------------------
    if (ConfPos["PLOT_HOR_STANDFORD"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["HPE"], PosIdx["HPL"], PosIdx["SOL"]])

        print('Plot Horizontal Stanford Diagram in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotHorStand(PosFile, PosData, Sol)

    # VERTICAL Stanford DIAGRAM
    # ----------------------------------------------------------
    if (ConfPos["PLOT_VER_STANDFORD"] == 1):
        # Read the cols we need from PosFile file
        PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None, \
                           usecols=[PosIdx["VPE"], PosIdx["VPL"], PosIdx["SOL"]])

        print('Plot Vertical Stanford Diagram in ' + Sol[0] + ' mode...')

        # Configure plot and call plot generation function
        plotVerStand(PosFile, PosData, Sol)