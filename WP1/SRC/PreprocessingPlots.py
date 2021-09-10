#!/usr/bin/env python

########################################################################
# PETRUS/SRC/PreprocessingPlots.py:
# This is the PreprocessingPlots Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           PreprocessingPlots.py
#  Date(YY/MM/DD): 05/02/21
#
#   Author: GNSS Academy
#   Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import PreproIdx
from InputOutput import REJECTION_CAUSE_DESC
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
from COMMON.Plots import generatePlot, generateChallengePlot
import numpy as np
from collections import OrderedDict

def initPlot(PreproObsFile, PlotConf, Title, Label):
    dirName = os.path.dirname(PreproObsFile)
    PreproObsFileName = os.path.basename(PreproObsFile)
    PreproObsFileNameSplit = PreproObsFileName.split('_')
    Rcvr = PreproObsFileNameSplit[2]
    DatepDat = PreproObsFileNameSplit[3]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PPVE/FIGURES/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)

# Plot Satellite Visibility
def plot_sat_visibility(PreproObsFile, PreproObsData):
    PlotConf = {}
    initPlot(PreproObsFile, PlotConf, "Satellite Visibility", "SAT_VISIBILITY")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4, 6.6)

    PlotConf["zLabel"] = "Elevation [deg]"

    PlotConf["yLabel"] = "GPS-PRN"
    PlotConf["yTicks"] = sorted(unique(PreproObsData[PreproIdx["PRN"]]))
    PlotConf["yTicksLabels"] = sorted(unique(PreproObsData[PreproIdx["PRN"]]))
    PlotConf["yLim"] = [0, max(unique(PreproObsData[PreproIdx["PRN"]]))]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False
    PlotConf["NotConv"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["MarkerSize"] = 5
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["xDataNotConv"] = {}
    PlotConf["yDataNotConv"] = {}
    PlotConf["zDataNotConv"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PRN"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 0
    PlotConf["xDataNotConv"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yDataNotConv"][Label] = PreproObsData[PreproIdx["PRN"]][FilterCond]
    PlotConf["zDataNotConv"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
# Plot Number of Satellites
def plot_num_sats(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Number of Satellites", "SAT_NUM")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4, 6.6)

    PlotConf["yLabel"] = "Number of Satellites"
    PlotConf["yTicks"] = range(0, 15, 2)
    PlotConf["yLim"] = [0, 14]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = True
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = '-'
    PlotConf["MarkerSize"] = 1
    PlotConf["LineWidth"] = 1

    RawSats = []
    SmoothedSats = []

    for sod in unique(PreproObsData[PreproIdx["SOD"]]):
        FilterCond = (PreproObsData[PreproIdx["SOD"]] == sod)
        RawSats.append(len(PreproObsData[PreproIdx["SOD"]][FilterCond]))

        FilterCond = ((PreproObsData[PreproIdx["SOD"]] == sod) & (PreproObsData[PreproIdx["STATUS"]] == 1))
        SmoothedSats.append(len(PreproObsData[PreproIdx["SOD"]][FilterCond]))

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Label"] = {}
    Label = 0
    PlotConf["Label"][Label] = 'Raw'
    PlotConf["Color"][Label] = 'orange'
    PlotConf["xData"][Label] = unique(PreproObsData[PreproIdx["SOD"]]) / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = RawSats
    Label = 1
    PlotConf["Label"][Label] = 'Smoothed'
    PlotConf["Color"][Label] = 'green'
    PlotConf["xData"][Label] = unique(PreproObsData[PreproIdx["SOD"]]) / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SmoothedSats

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
# Plot Satellite Polar View
def plot_sat_polar_view(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Satellite Polar View", "SAT_POLAR_VIEW")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10,10)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32

    PlotConf["tData"] = {}
    PlotConf["rData"] = {}
    PlotConf["zData"] = {}
    PlotConf["tData"] = np.deg2rad(PreproObsData[PreproIdx["AZIM"]])
    PlotConf["rData"] = PreproObsData[PreproIdx["ELEV"]]
    PlotConf["zData"] = PreproObsData[PreproIdx["PRN"]]

    # Call generatePlot from Plots library
    generateChallengePlot(PlotConf)

# Plot C1 - C1Smoothed
def plot_C1_C1smooth_diff(PreproObsFile, PreproObsData):

    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "C1 - C1Smoothed", "SAT_C1_C1SMOOTHED")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["yLabel"] = "C1 - C1Smoothed [m]"
    PlotConf["yTicks"] = range(-3, 3)
    PlotConf["yLim"] = [-3, 2.25]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "C/N0 [deg]"
    PlotConf["ColorBarMin"] = 35.
    PlotConf["ColorBarMax"] = 51.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["C1"]][FilterCond] - PreproObsData[PreproIdx["C1SMOOTHED"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["S1"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
def plot_C1_C1smooth_diff_vs_elev(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "C1 - C1Smoothed vs Elevation", "SAT_C1_C1SMOOTHED_VS_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["yLabel"] = "C1 - C1Smoothed [m]"
    PlotConf["yTicks"] = range(-3, 3)
    PlotConf["yLim"] = [-3, 2.25]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "C/N0 [deg]"
    PlotConf["ColorBarMin"] = 35.
    PlotConf["ColorBarMax"] = 51.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["C1"]][FilterCond] - PreproObsData[PreproIdx["C1SMOOTHED"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["S1"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Rejection Flags
def plot_rejection_flags(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Rejection Flags", "REJ_FLAGS")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "Rejection Flags"
    PlotConf["yTicks"] = range(1, 11)
    PlotConf["yTicksLabels"] = REJECTION_CAUSE_DESC.keys()
    PlotConf["yLim"] = [0, 11]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False
    PlotConf["RejectFlag"] = True

    PlotConf["ColorBar"] = "gist_ncar"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 32.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["REJECT"]] != 0
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["REJECT"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Code Rates
def plot_code_rate(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Code Rate", "CODE_RATE")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "Code Rate [m/s]"
    PlotConf["yTicks"] = range(-800, 1000, 200)
    PlotConf["yLim"] = [-800, 800]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["CODE RATE"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
def plot_code_acc(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Code Rate Step", "CODE_RATE_STEP")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "Code Rate Step [m/s^2]"
    PlotConf["yTicks"] = np.arange(-0.05, 0.25, 0.05)
    PlotConf["yLim"] = [-0.05, 0.20]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["CODE ACC"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Phase Rates
def plot_phase_rate(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Phase Rate", "PHASE_RATE")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "Phase Rate [m/s]"
    PlotConf["yTicks"] = range(-800, 1000, 200)
    PlotConf["yLim"] = [-800, 800]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PHASE RATE"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
def plot_phase_acc(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Phase Rate Step", "PHASE_RATE_STEP")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "Phase Rate Step [m/s^2]"
    PlotConf["yTicks"] = np.arange(-0.05, 0.25, 0.05)
    PlotConf["yLim"] = [-0.05, 0.20]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PHASE ACC"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# VTEC Gradient
def plot_VTEC_gradient(PreproObsFile, PreproObsData):
    PlotConf = {}
    initPlot(PreproObsFile, PlotConf, "VTEC Gradient", "VTEC_RATE")
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "VTEC Gradient [mm/s]"
    PlotConf["yTicks"] = range(-4, 5)
    PlotConf["yLim"] = [-4, 4]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond]  / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["VTEC RATE"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
# AATR index
def plot_AATR(PreproObsFile, PreproObsData):

    PlotConf = {}
    initPlot(PreproObsFile, PlotConf, "AATR", "AATR")
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)

    PlotConf["yLabel"] = "AATR [mm/s]"
    PlotConf["yTicks"] = range(-4, 5)
    PlotConf["yLim"] = [-4, 4]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = 'P'
    PlotConf["MarkerSize"] = 0.5
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    FilterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond]  / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["iAATR"]][FilterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def generatePreproPlots(PreproObsFile):
    # Purpose: generate output plots regarding Preprocessing results

    # Parameters
    # ==========
    # PreproObsFile: str
    #         Path to PREPRO OBS output file

    # Returns
    # =======
    # Nothing

    # Satellite Visibility
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,
    usecols=[PreproIdx["SOD"], PreproIdx["PRN"], PreproIdx["ELEV"], PreproIdx["STATUS"]])

    plot_sat_visibility(PreproObsFile, PreproObsData)

    # Number of Satellites
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"]])

    plot_num_sats(PreproObsFile, PreproObsData)

    # Satellite Polar View
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["PRN"], PreproIdx["ELEV"], PreproIdx["AZIM"]])

    plot_sat_polar_view(PreproObsFile, PreproObsData)

    # Satellite C1 - C1Smoothed
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["C1"], PreproIdx["C1SMOOTHED"], PreproIdx["S1"]])

    plot_C1_C1smooth_diff(PreproObsFile, PreproObsData)

    # Satellite C1 - C1Smoothed vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["C1"], PreproIdx["C1SMOOTHED"], PreproIdx["S1"]])

    plot_C1_C1smooth_diff_vs_elev(PreproObsFile, PreproObsData)

    # Satellite Rejection Flags
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["PRN"], PreproIdx["REJECT"]])

    plot_rejection_flags(PreproObsFile, PreproObsData)

    # Satellite Code Rate
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["CODE RATE"], PreproIdx["ELEV"]])

    plot_code_rate(PreproObsFile, PreproObsData)

    # Satellite Code Rate Step
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["CODE ACC"], PreproIdx["ELEV"]])

    plot_code_acc(PreproObsFile, PreproObsData)

    # Satellite Phase Rate
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["PHASE RATE"], PreproIdx["ELEV"]])

    plot_phase_rate(PreproObsFile, PreproObsData)

    # Satellite Phase Rate Step
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["PHASE ACC"], PreproIdx["ELEV"]])

    plot_phase_acc(PreproObsFile, PreproObsData)

    # Satellite VTEC Gradient
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["VTEC RATE"], PreproIdx["ELEV"]])

    plot_VTEC_gradient(PreproObsFile, PreproObsData)

    # Satellite Instantaneus AATR
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"], PreproIdx["iAATR"], PreproIdx["ELEV"]])

    plot_AATR(PreproObsFile, PreproObsData)

if __name__ == "__main__":
    generatePreproPlots('/Users/jorialand/Desktop/PROJECTS/GNSS-ACADEMY/PETRUS/WP1/SCN/SCEN-PETRUS-GEO123-JAN19/OUT/PPVE/PREPRO_OBS_TLSZ_Y19D014.dat')
    generatePreproPlots('/Users/jorialand/Desktop/PROJECTS/GNSS-ACADEMY/PETRUS/WP1/SCN/SCEN-PETRUS-GEO123-JAN19/OUT/PPVE/PREPRO_OBS_TLSA_Y19D014.dat')
    generatePreproPlots('/Users/jorialand/Desktop/PROJECTS/GNSS-ACADEMY/PETRUS/WP1/SCN/SCEN-PETRUS-GEO123-JAN19/OUT/PPVE/PREPRO_OBS_RKKA_Y19D014.dat')
    generatePreproPlots('/Users/jorialand/Desktop/PROJECTS/GNSS-ACADEMY/PETRUS/WP1/SCN/SCEN-PETRUS-GEO123-JAN19/OUT/PPVE/PREPRO_OBS_LPIA_Y19D014.dat')
    generatePreproPlots('/Users/jorialand/Desktop/PROJECTS/GNSS-ACADEMY/PETRUS/WP1/SCN/SCEN-PETRUS-GEO123-JAN19/OUT/PPVE/PREPRO_OBS_LAPA_Y19D014.dat')
    generatePreproPlots('/Users/jorialand/Desktop/PROJECTS/GNSS-ACADEMY/PETRUS/WP1/SCN/SCEN-PETRUS-GEO123-JAN19/OUT/PPVE/PREPRO_OBS_ALYA_Y19D014.dat')