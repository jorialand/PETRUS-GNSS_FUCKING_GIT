
import sys, os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import conda
CondaFileDir = conda.__file__
CondaDir = CondaFileDir.split('lib')[0]
ProjLib = os.path.join(os.path.join(CondaDir, 'share'), 'proj')
os.environ["PROJ_LIB"] = ProjLib
from mpl_toolkits.basemap import Basemap
from scipy.stats import norm

import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

import PlotsConstants as Const

def createFigure(PlotConf, Polar = 0):
    if Polar == 0:
        if "FigSize" in PlotConf:
            fig, ax = plt.subplots(1, 1, figsize = PlotConf["FigSize"])
        
        else:
            fig, ax = plt.subplots(1, 1)

    else:
        if "FigSize" in PlotConf:
            fig = plt.figure(figsize = PlotConf["FigSize"])
            ax = fig.add_subplot(111, polar=True)

        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, polar=True)

    return fig, ax

def saveFigure(fig, Path):
    Dir = os.path.dirname(Path)
    try:
        os.makedirs(Dir)
    except: pass

    fig.savefig(Path, dpi=200., bbox_inches='tight')

def prepareAxis(PlotConf, ax):
    for key in PlotConf:
        if key == "Title":
            ax.set_title(PlotConf["Title"])

        for axis in ["x", "y"]:
            if axis == "x":
                if key == axis + "Label":
                    ax.set_xlabel(PlotConf[axis + "Label"])

                if key == axis + "Ticks":
                    ax.set_xticks(PlotConf[axis + "Ticks"])

                if key == axis + "TicksLabels":
                    ax.set_xticklabels(PlotConf[axis + "TicksLabels"])
                
                if key == axis + "Lim":
                    ax.set_xlim(PlotConf[axis + "Lim"])

            if axis == "y":
                if key == axis + "Label":
                    ax.set_ylabel(PlotConf[axis + "Label"])

                if key == axis + "Ticks":
                    ax.set_yticks(PlotConf[axis + "Ticks"])

                if key == axis + "TicksLabels":
                    ax.set_yticklabels(PlotConf[axis + "TicksLabels"])
                
                if key == axis + "Lim":
                    ax.set_ylim(PlotConf[axis + "Lim"])

        if key == "Grid" and PlotConf[key] == True:
            ax.grid(linestyle='--', linewidth=0.5, which='both')

    try:
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_scientific(False)
    except:
        pass

def prepareColorBar(PlotConf, ax, Values, Polar = 0):
    try:
        Min = PlotConf["ColorBarMin"]
    except:
        Mins = []
        for v in Values.values():
            Mins.append(min(v))
        Min = min(Mins)
    try:
        Max = PlotConf["ColorBarMax"]
    except:
        Maxs = []
        for v in Values.values():
            Maxs.append(max(v))
        Max = max(Maxs)
    normalize = mpl.cm.colors.Normalize(vmin=Min, vmax=Max)

    if Polar==0:
        divider = make_axes_locatable(ax)

        color_ax = divider.append_axes("right", size="3%", pad="2%")

    else:
        color_ax, kwargs = mpl.colorbar.make_axes(ax, orientation="vertical", pad=0.05, fraction=0.05)
    
    cmap = mpl.cm.get_cmap(PlotConf["ColorBar"])
    
    if "ColorBarTicks" in PlotConf:
        cbar = mpl.colorbar.ColorbarBase(color_ax, 
        cmap=cmap,
        norm=mpl.colors.BoundaryNorm(PlotConf["ColorBarTicks"], cmap.N),
        label=PlotConf["ColorBarLabel"],
        ticks=PlotConf["ColorBarTicks"])

    else:
        cbar = mpl.colorbar.ColorbarBase(color_ax, 
        cmap=cmap,
        norm=mpl.colors.Normalize(vmin=Min, vmax=Max),
        label=PlotConf["ColorBarLabel"])

    return normalize, cmap

def drawMap(PlotConf, ax,):
    Map = Basemap(projection = 'cyl',
    llcrnrlat  = PlotConf["LatMin"]-0,
    urcrnrlat  = PlotConf["LatMax"]+0,
    llcrnrlon  = PlotConf["LonMin"]-0,
    urcrnrlon  = PlotConf["LonMax"]+0,
    lat_ts     = 10,
    resolution = 'l',
    ax         = ax)
    
    if not ("xTicks" in PlotConf or "yTicks" in PlotConf):
        # Draw map meridians
        Map.drawmeridians(
        np.arange(PlotConf["LonMin"],PlotConf["LonMax"]+1,PlotConf["LonStep"]),
        labels = [0,0,0,1],
        fontsize = 6,
        linewidth=0.2)
            
        # Draw map parallels
        Map.drawparallels(
        np.arange(PlotConf["LatMin"],PlotConf["LatMax"]+1,PlotConf["LatStep"]),
        labels = [1,0,0,0],
        fontsize = 6,
        linewidth=0.2)

    # Draw coastlines
    Map.drawcoastlines(linewidth=0.5)

    # Draw countries
    Map.drawcountries(linewidth=0.25)

def generateLinesPlot(PlotConf):
    fig, ax = createFigure(PlotConf)

    for key in PlotConf:
        if key == "ColorBar":
            normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"])
        if key == "Map" and PlotConf[key] == True:
            drawMap(PlotConf, ax)

    for Label in PlotConf["yData"].keys():
        LineWidth = 1.5
        Alpha = 1

        for key in PlotConf:
            if key == "LineWidth":
                LineWidth = PlotConf["LineWidth"]
            if key == "Alpha":
                Alpha = PlotConf["Alpha"]

        if "MarkerSize" in PlotConf and Label in PlotConf["MarkerSize"]:
            LineWidth = PlotConf["MarkerSize"][Label]

        if "ColorBar" in PlotConf:
            if "Color" in PlotConf and Label in PlotConf["Color"]:
                ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label],
                marker = PlotConf["Marker"],
                linewidth = LineWidth,
                s = LineWidth,
                alpha=Alpha,
                c = PlotConf["Color"][Label])

            else:
                ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label],
                marker = PlotConf["Marker"],
                linewidth = LineWidth,
                s = LineWidth,
                alpha=Alpha,
                c = cmap(normalize(np.array(PlotConf["zData"][Label]))))

        else:
            ax.plot(PlotConf["xData"][Label], PlotConf["yData"][Label],
            PlotConf["Marker"],
            linewidth = LineWidth,
            markersize = LineWidth,
            alpha=Alpha,
            c = PlotConf["Color"][Label])

    for key in PlotConf:
        if key == "Legend" and PlotConf[key] == 1:
            Legend = plt.legend(PlotConf["yData"].keys())
            for iHandle in range(len(Legend.legendHandles)):
                Legend.legendHandles[iHandle]._legmarker.set_markersize(7)

    prepareAxis(PlotConf, ax)

    saveFigure(fig, PlotConf["Path"])
    plt.close('all')

def generatePolarPlot(PlotConf):
    LineWidth = 1.5

    fig, ax = createFigure(PlotConf, 1)

    for key in PlotConf:

        if key == "Title":
            ax.set_title(PlotConf["Title"])

        if key == "Grid" and PlotConf[key] == True:
            ax.grid(linestyle='--', linewidth=0.5, which='both')

        if key == "LineWidth":
            LineWidth = PlotConf["LineWidth"]

        if key == "ColorBar":
            normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"], 1)

    for Label in PlotConf["rData"].keys():
        if "ColorBar" in PlotConf:
            ax.scatter(PlotConf["tData"][Label], PlotConf["rData"][Label],
            marker = PlotConf["Marker"],
            linewidth = LineWidth,
            s = LineWidth,
            c = cmap(normalize(np.array(PlotConf["zData"][Label]))))

        else:
            ax.plot(PlotConf["tData"][Label], PlotConf["rData"][Label], 'ro',
            marker = PlotConf["Marker"],
            linewidth = LineWidth,
            markersize = LineWidth)

        ax.set_rmin(90)
        ax.set_rmax(0)
        ax.set_xticks([0, np.pi/2, np.pi, (3/2)*np.pi])
        ax.set_xticklabels(['E', 'N', 'W', 'S'])

    saveFigure(fig, PlotConf["Path"])
    plt.close('all')

def generateStatsPlot(PlotConf):
    BarWidth = 0.3

    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)

    CellText = []
    StatsList = list(list(PlotConf["yData"].values())[0].keys())
    NStats = len(StatsList)
    RowColors = plt.cm.rainbow(np.linspace(0.3, 0.9, NStats))[::-1]
    ColColors = plt.cm.binary(np.linspace(0.2, 0.2, len(PlotConf["xData"])))
    Colors = {}

    for iLabel, Label in enumerate(PlotConf["xData"]):
        Colors[iLabel] = {}

        for iStat, Stat in enumerate(PlotConf["yData"][Label].values()):
            Colors[iLabel][iStat] = iStat

            ax.bar(iLabel, Stat, 
            width = BarWidth,
            color=RowColors[Colors[iLabel][iStat]])

    for Stat in StatsList:
        List = []
        for iLabel, Label in enumerate(PlotConf["xData"]):
            List.append("%.2f" % PlotConf["yData"][Label][Stat])
        CellText.append(List)

    ax.set_xticks(range(len(PlotConf["xData"])))
    for tick in ax.xaxis.get_major_ticks():
        # tick.tick1On = tick.tick2On = False
        # tick.label1On = tick.label2On = False	
        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False) 
        tick.label1.set_visible(False)
        tick.label2.set_visible(False)

    ax.set_xlim(left=-.5, right=len(PlotConf["xData"]) - 1 + .5)

    ax.table(cellText=CellText,
    rowLabels=StatsList,
    rowColours=RowColors,
    colLabels=PlotConf["xData"],
    colColours=ColColors,
    cellLoc='center',
    loc='bottom',
    fontsize='small')

    plt.subplots_adjust(left=0.22, bottom=0.2)

    saveFigure(fig, PlotConf["Path"])
    plt.close('all')

def generatePlot(PlotConf):
    if(PlotConf["Type"] == "Lines"):
        generateLinesPlot(PlotConf)
    
    if(PlotConf["Type"] == "Polar"):
        generatePolarPlot(PlotConf)
    
    if(PlotConf["Type"] == "Stats"):
        generateStatsPlot(PlotConf)
