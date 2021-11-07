
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

import COMMON.PlotsConstants as Const


def saveFigure(fig, Path):
    Dir = os.path.dirname(Path)
    try:
        os.makedirs(Dir)
    except:
        pass
    fig.savefig(Path, dpi=150., bbox_inches='tight')


def createFigure(PlotConf):
    # Delete previous opened figures, to avoid performance issues
    plt.close('all')
    try:
        fig, ax = plt.subplots(1, 1, figsize=PlotConf["FigSize"])

    except:
        fig, ax = plt.subplots(1, 1)

    return fig, ax


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


def prepareColorBar(PlotConf, ax, Values):
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

    divider = make_axes_locatable(ax)
    # size of the plot and gap of pad from the plot
    color_ax = divider.append_axes("right", size="3%", pad="2%")
    cmap = mpl.cm.get_cmap(PlotConf["ColorBar"])
    cbar = mpl.colorbar.ColorbarBase(color_ax,
                                     cmap=cmap,
                                     norm=mpl.colors.Normalize(vmin=Min, vmax=Max),
                                     ticks=PlotConf["ColorBarTicks"],
                                     label=PlotConf["ColorBarLabel"])

    return normalize, cmap


def drawMap(PlotConf, ax, ):
    Map = Basemap(projection='cyl',
                  llcrnrlat=PlotConf["LatMin"] - 0,
                  urcrnrlat=PlotConf["LatMax"] + 0,
                  llcrnrlon=PlotConf["LonMin"] - 0,
                  urcrnrlon=PlotConf["LonMax"] + 0,
                  lat_ts=10,
                  resolution='l',
                  ax=ax)

    # Draw map meridians
    Map.drawmeridians(
        np.arange(PlotConf["LonMin"], PlotConf["LonMax"] + 1, PlotConf["LonStep"]),
        labels=[0, 0, 0, 1],
        fontsize=6,
        linewidth=0.2)

    # Draw map parallels
    Map.drawparallels(
        np.arange(PlotConf["LatMin"], PlotConf["LatMax"] + 1, PlotConf["LatStep"]),
        labels=[1, 0, 0, 0],
        fontsize=6,
        linewidth=0.2)

    # Draw coastlines
    Map.drawcoastlines(linewidth=0.5)

    # Draw countries
    Map.drawcountries(linewidth=0.25)


def generateLinesPlot(PlotConf):
    LineWidth = 1.5

    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)

    for key in PlotConf:
        if key == "LineWidth":
            LineWidth = PlotConf["LineWidth"]
        if key == "ColorBar":
            normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"])
        if key == "Map" and PlotConf[key] == True:
            drawMap(PlotConf, ax)

    for Label in PlotConf["yData"].keys():
        if "Background" in PlotConf and Label in PlotConf["Background"]:
            ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label],
                       marker=PlotConf["BackMarker"],
                       s=PlotConf["BackLineWidth"],
                       zorder=1,
                       c=PlotConf["ColorMarker"] if "ColorMarker" in PlotConf else None)
        elif "SecondAxis" in PlotConf and Label in PlotConf["SecondAxis"]:
            ax2 = ax.twinx()
            ax2.set_ylabel(PlotConf["y2Label"])
            ax2.set_ylim(PlotConf["y2Lim"][0], PlotConf["y2Lim"][1])
            ax2.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label],
                        marker=PlotConf["Marker"],
                        s=LineWidth,
                        c="black",
                        label=Label,
                        zorder=10)
            ln2, lab2 = ax2.get_legend_handles_labels()
        elif "ColorBar" in PlotConf:
            ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label],
                       marker=PlotConf["Marker"],
                       s=LineWidth,
                       zorder=10,
                       c=cmap(normalize(np.array(PlotConf["zData"][Label]))))
        else:
            ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label],
                       marker=PlotConf["Marker"],
                       s=LineWidth,
                       label=Label,
                       zorder=10)
            ln1, lab1 = ax.get_legend_handles_labels()

    if "HLine" in PlotConf.keys():
        for Line in PlotConf["HLine"]:
            ax.hlines(Line[0], Line[1], Line[2], color="black", linestyle='--', linewidth=0.75)
    if "VLine" in PlotConf.keys():
        for Line in PlotConf["VLine"]:
            ax.vlines(Line[0], Line[1], Line[2], color="black", linestyle='--', linewidth=0.75)
    if "SLine" in PlotConf.keys():
        ax.plot(PlotConf["SLine"][0], PlotConf["SLine"][1], color="black", linestyle='--', linewidth=0.75)

    if "Legend" in PlotConf.keys():
        if "SecondAxis" not in PlotConf:
            plt.legend(PlotConf["Legend"], markerscale=4.0, loc="upper right")
        elif "SecondAxis" in PlotConf:
            ln = ln1 + ln2
            lab = lab1 + lab2
            ax.legend(ln, lab, markerscale=4.0, loc="upper right")

    saveFigure(fig, PlotConf["Path"])


def generateBarsPlot(PlotConf):
    BarWidth = 1.0

    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)
    ax.set_xticks([])

    for key in PlotConf:
        if key == "BarWidth":
            BarWidth = PlotConf["BarWidth"]

    TableInfo = []
    for Label in PlotConf["yData"].keys():
        ax.bar(PlotConf["xData"][Label], PlotConf["yData"][Label], width=BarWidth)
        TableInfo.append(PlotConf["yData"][Label])

    Rows = PlotConf["Legend"]
    Cols = ["G" + "%02d" % i for i in range(1, 33, 1)]
    ax.table(cellText=TableInfo, rowLabels=Rows, colLabels=Cols, loc="bottom")

    for key in PlotConf:
        if key == "Legend":
            plt.legend(PlotConf["Legend"])

    saveFigure(fig, PlotConf["Path"])


def generateHistogram(PlotConf):
    BarWidth = 0.005

    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)

    for key in PlotConf:
        if key == "BarWidth":
            BarWidth = PlotConf["BarWidth"]

    ax.bar(PlotConf["xData"], PlotConf["yData"], width=BarWidth, align="edge")
    ax.plot(PlotConf["x2Data"], PlotConf["y2Data"], c="orange")

    for key in PlotConf:
        if key == "Legend":
            plt.legend(PlotConf["Legend"])

    saveFigure(fig, PlotConf["Path"])


def generateMapPlot(PlotConf):
    LineWidth = 15.0

    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)

    for key in PlotConf:
        if key == "LineWidth":
            LineWidth = PlotConf["LineWidth"]
        if key == "ColorBar":
            normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"])
        if key == "Map" and PlotConf[key] == True:
            drawMap(PlotConf, ax)

    ax.scatter(PlotConf["xData"], PlotConf["yData"],
               marker=PlotConf["Marker"],
               s=LineWidth,
               zorder=10,
               c=cmap(normalize(np.array(PlotConf["zData"]))))

    for i, label1 in enumerate(PlotConf["nData"]):
        ax.annotate("%s" % label1, xy=(PlotConf["xData"][i], PlotConf["yData"][i]),
                    xytext=(PlotConf["xData"][i] - 2.0, (PlotConf["yData"][i] + 0.7)))

    for j, label2 in enumerate(PlotConf["zData"]):
        if "Decimal" in PlotConf.keys():
            ax.annotate("%.1e" % label2, xy=(PlotConf["xData"][j], PlotConf["yData"][j]),
                        xytext=(PlotConf["xData"][j] - 2.0, (PlotConf["yData"][j] - 1.5)))
        elif "Integer" in PlotConf.keys():
            ax.annotate("%d" % label2, xy=(PlotConf["xData"][j], PlotConf["yData"][j]),
                        xytext=(PlotConf["xData"][j] - 1.0, (PlotConf["yData"][j] - 1.5)))
        else:
            ax.annotate("%.2f" % label2, xy=(PlotConf["xData"][j], PlotConf["yData"][j]),
                        xytext=(PlotConf["xData"][j] - 2.0, (PlotConf["yData"][j] - 1.5)))

    saveFigure(fig, PlotConf["Path"])


def generatePlot(PlotConf):
    if (PlotConf["Type"] == "Lines"):
        generateLinesPlot(PlotConf)
    elif (PlotConf["Type"] == "Bars"):
        generateBarsPlot(PlotConf)
    elif (PlotConf["Type"] == "Hist"):
        generateHistogram(PlotConf)
    elif (PlotConf["Type"] == "Map"):
        generateMapPlot(PlotConf)
