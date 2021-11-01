#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Perf.py:
# This is the Performances Module of PETRUS tool
#
# Project:        PETRUS
# File:           Perf.py
# Date(YY/MM/DD): 16/02/21
#
# Author: GNSS Academy
# Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

# Import External and Internal functions and Libraries
# ----------------------------------------------------------------------
import sys, os

from scipy import stats

# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from InputOutput import RcvrIdx
from COMMON import Stats, GnssConstants
from math import sqrt
from collections import OrderedDict
from InputOutput import generateHistFile


def UpdateBuff(ContBuff, Status, Epochs):
    for Epoch in range(int(Epochs)):
        ContBuff.pop(0)
        ContBuff.append(Status)

def initializePerfInfo(Conf, Services, Rcvr, RcvrInfo, Doy, PerfInfo, VpeHistInfo):
    """
    Initialize the dict for the performances information.
    Initialize the dict for the LPV200 histogram.
    :param Conf:
    :param Services:
    :param Rcvr:
    :param RcvrInfo:
    :param Doy:
    :param PerfInfo:
    :param VpeHistInfo:
    :return:
    """

    # Initialize internal variables
    Idx = {"FLAG": 0, "HAL": 1, "VAL": 2, "HPE95": 3, "VPE95": 4, "VPE1E7": 5, "AVAI": 6, "CONT": 7, "CINT": 8}

    # Loop over all the activated service levels
    for Service in Services:
        # If service activated
        if int(Conf[Service][Idx["FLAG"]]) == 1:
            # Initialize PerInfo dictionary
            PerfInfo[Service] = {
                "Rcvr": Rcvr,  # Receiver acronym
                "Lon": float(RcvrInfo[RcvrIdx["LON"]]),  # Receiver reference longitude
                "Lat": float(RcvrInfo[RcvrIdx["LAT"]]),  # Receiver reference latitude
                "Doy": Doy,  # Day of year
                "Service": Service,  # Service level
                "SamSol": 86400 // int(Conf["SAMPLING_RATE"]),  # Number of total samples processed
                "SamNoSol": 86400 // int(Conf["SAMPLING_RATE"]),  # Number of samples with no SBAS solution
                "Avail": 0,  # Availability percentage
                "ContRisk": 0.0,  # Continuity risk
                "ContBuff": [0] * int(Conf[Service][Idx["CINT"]]),  # Continuity risk buffer
                "PrevStatus": 0,  # Previous availability status
                "PrevSod": 0.0,  # Previous computed epoch
                "ContEvent": 0,  # Number of discontinuity events
                "NotAvail": 0,  # Number of non-available samples for the selected service level
                "NsvMin": 1000,  # Minimum number of satellites
                "NsvMax": 0,  # Maximum number of satellites
                "HpeRms": 0.0,  # HPE RMS
                "VpeRms": 0.0,  # VPE RMS
                "Hpe95": 0.0,  # HPE 95% percentile
                "HpeHist": {},  # HPE Histogram
                "Vpe95": 0.0,  # VPE 95% percentile
                "VpeHist": {},  # VPE Histogram
                "HpeMax": 0.0,  # Maximum HPE
                "VpeMax": 0.0,  # Maximum VPE
                "ExtVpe": 0.0,  # Extrapolated VPE
                "HplMin": 1000.0,  # Minimum HPL
                "VplMin": 1000.0,  # Minimum VPL
                "HplMax": 0.0,  # Maximum HPL
                "VplMax": 0.0,  # Maximum VPL
                "HsiMax": 0.0,  # Maximum HSI
                "VsiMax": 0.0,  # Maximum VSI
                "Nmi": 0,  # Number of misleading information events
                "Nhmi": 0,  # Number of hazardous misleading information events
                "PdopMax": 0.0,  # Maximum PDOP
                "HdopMax": 0.0,  # Maximum HDOP
                "VdopMax": 0.0,  # Maximum VDOP
            }  # End of PerfInfo[Service]

            if Service == "LPV200":
                # Initialize LPV200 VpeHistInfo dictionary
                VpeHistInfo["Rcvr"] = Rcvr  # Receiver acronym
                VpeHistInfo["Service"] = Service  # Service level
                VpeHistInfo["BinId"] = 0  # Bin ID
                VpeHistInfo["BinMin"] = 0.0  # Bin minimum
                VpeHistInfo["BinMax"] = 0.0  # Bin maximum
                VpeHistInfo["BinNumSam"] = 0  # Bin number of samples
                VpeHistInfo["BinFreq"] = 0.0  # Bin relative frequency

    # Check if there are active service levels
    if len(PerfInfo.keys()) == 0:
        sys.stderr.write("ERROR: Please activate at least one service level in the configuration file \n")
        sys.exit(1)

def UpdatePerfEpoch(Conf, Service, PosInfo, PerfInfoSer):
    """
    Update the performance info, this includes Protection Levels, Safety Indexes, DOPs, Availability, Position Errors,
    and Continuity
    :param Conf:
    :param Service:
    :param PosInfo:
    :param PerfInfoSer:
    :return:
    """

    # Initialize internal variables
    HistRes = GnssConstants.HIST_RES
    AvailStatus = 0
    Idx = {"FLAG": 0, "HAL": 1, "VAL": 2, "HPE95": 3, "VPE95": 4, "VPE1E7": 5, "AVAI": 6, "CONT": 7, "CINT": 8}

    # If SBAS solution has been achieved
    if PosInfo["Sol"] != 0:

        # Measurements with no SBAS solution
        PerfInfoSer["SamNoSol"] = PerfInfoSer["SamNoSol"] - 1

        # SVs used
        PerfInfoSer["NsvMin"] = Stats.UpdateMin(PerfInfoSer["NsvMin"], PosInfo["NumSatSol"])
        PerfInfoSer["NsvMax"] = Stats.UpdateMax(PerfInfoSer["NsvMax"], PosInfo["NumSatSol"])

        # Protection levels
        PerfInfoSer["HplMin"] = Stats.UpdateMin(PerfInfoSer["HplMin"], PosInfo["Hpl"])
        PerfInfoSer["VplMin"] = Stats.UpdateMin(PerfInfoSer["VplMin"], PosInfo["Vpl"])
        PerfInfoSer["HplMax"] = Stats.UpdateMax(PerfInfoSer["HplMax"], PosInfo["Hpl"])
        PerfInfoSer["VplMax"] = Stats.UpdateMax(PerfInfoSer["VplMax"], PosInfo["Vpl"])

        # Safety Indexes
        PerfInfoSer["HsiMax"] = Stats.UpdateMax(PerfInfoSer["HsiMax"], abs(PosInfo["Hsi"]))
        PerfInfoSer["VsiMax"] = Stats.UpdateMax(PerfInfoSer["VsiMax"], abs(PosInfo["Vsi"]))

        # DOPs
        PerfInfoSer["PdopMax"] = Stats.UpdateMax(PerfInfoSer["PdopMax"], PosInfo["Pdop"])
        PerfInfoSer["HdopMax"] = Stats.UpdateMax(PerfInfoSer["HdopMax"], PosInfo["Hdop"])
        PerfInfoSer["VdopMax"] = Stats.UpdateMax(PerfInfoSer["VdopMax"], PosInfo["Vdop"])

        # Availability
        # PL > AL --> NOT-AVAILABLE
        if (PosInfo["Hpl"] / Conf[Service][Idx["HAL"]]) > 1 or (PosInfo["Vpl"] / Conf[Service][Idx["VAL"]]) > 1:
            PerfInfoSer["NotAvail"] = PerfInfoSer["NotAvail"] + 1
        # PL <= AL --> AVAILABLE
        elif (PosInfo["Hpl"] / Conf[Service][Idx["HAL"]]) <= 1 and (PosInfo["Vpl"] / Conf[Service][Idx["VAL"]]) <= 1:
            AvailStatus = 1
            PerfInfoSer["Avail"] = PerfInfoSer["Avail"] + 1

            # PE RMS
            PerfInfoSer["HpeRms"] = PerfInfoSer["HpeRms"] + PosInfo["Hpe"] ** 2
            PerfInfoSer["VpeRms"] = PerfInfoSer["VpeRms"] + PosInfo["Vpe"] ** 2

            # PE Histograms
            Stats.UpdateHist(PerfInfoSer["HpeHist"], abs(PosInfo["Hpe"]), HistRes)
            Stats.UpdateHist(PerfInfoSer["VpeHist"], abs(PosInfo["Vpe"]), HistRes)

            # PE
            PerfInfoSer["HpeMax"] = Stats.UpdateMax(PerfInfoSer["HpeMax"], PosInfo["Hpe"])
            PerfInfoSer["VpeMax"] = Stats.UpdateMax(PerfInfoSer["VpeMax"], PosInfo["Vpe"])

            # SI >= 1 --> MI - Misleading Information
            if PosInfo["Hsi"] >= 1 or abs(PosInfo["Vsi"]) >= 1:
                # MI - Misleading Information
                if PosInfo["Hpe"] < Conf[Service][Idx["HAL"]] and abs(PosInfo["Vpe"]) < Conf[Service][Idx["VAL"]]:
                    PerfInfoSer["Nmi"] = PerfInfoSer["Nmi"] + 1
                # HMI - Hazardous MI
                elif PosInfo["Hpe"] >= Conf[Service][Idx["HAL"]] or abs(PosInfo["Vpe"]) >= Conf[Service][Idx["VAL"]]:
                    PerfInfoSer["Nhmi"] = PerfInfoSer["Nhmi"] + 1


    # CR - Continuity Risk
    # Discontinuity Event --> n(AVAILABLE -> NOT-AVAILABLE)
    if AvailStatus == 0 and PerfInfoSer["PrevStatus"] == 1:
        PerfInfoSer["ContEvent"] = PerfInfoSer["ContEvent"] + sum(PerfInfoSer["ContBuff"])
    # DATA GAP Detected
    gap = PosInfo["Sod"] - PerfInfoSer["PrevSod"]
    if PerfInfoSer["PrevSod"] != 0. and gap > int(Conf["SAMPLING_RATE"]):
        PerfInfoSer["ContEvent"] = PerfInfoSer["ContEvent"] + sum(PerfInfoSer["ContBuff"])
        # Include gap in the continuity buffer if the Hatch Filter has not been reset
        if gap < int(Conf["HATCH_GAP_TH"]):
            UpdateBuff(PerfInfoSer["ContBuff"], 0, gap)
        # Reset the continuity buffer if the Hatch Filter has been reset
        elif gap >= int(Conf["HATCH_GAP_TH"]):
            PerfInfoSer["ContBuff"] = [0] * int(Conf[Service][Idx["CINT"]])

    # Update continuity buffer with current availability status
    UpdateBuff(PerfInfoSer["ContBuff"], AvailStatus, 1)
    # Update previous availabilty status and previous computed epoch
    PerfInfoSer["PrevStatus"] = AvailStatus
    PerfInfoSer["PrevSod"] = PosInfo["Sod"]

def computeFinalPerf(PerfInfoSer):
    """
    Compute the performances once all the observations are processed.

    :param PerfInfoSer:
    :return:
    """

    # Number of non-available samples for this Service Level
    PerfInfoSer["NotAvail"] = PerfInfoSer["NotAvail"] + PerfInfoSer["SamNoSol"]

    # At least a sample with availability
    if PerfInfoSer["Avail"] > 0:

        # Compute final HPE and VPE RMS
        # ----------------------------------------------------------------------
        # Compute HPE RMS
        PerfInfoSer["HpeRms"] = sqrt(PerfInfoSer["HpeRms"] / PerfInfoSer["Avail"])
        # Compute VPE RMS
        PerfInfoSer["VpeRms"] = sqrt(PerfInfoSer["VpeRms"] / PerfInfoSer["Avail"])

        # Compute final HPE95 and VPE95 values
        # ----------------------------------------------------------------------
        # HPE95
        Cdf, Sigmas = Stats.computeCdfFromHistogram(PerfInfoSer["HpeHist"], PerfInfoSer["Avail"])
        PerfInfoSer["Hpe95"] = Stats.computePercentile(Cdf, 95)
        # VPE95
        Cdf, Sigmas = Stats.computeCdfFromHistogram(PerfInfoSer["VpeHist"], PerfInfoSer["Avail"])
        PerfInfoSer["Vpe95"] = Stats.computePercentile(Cdf, 95)

        # Compute final extrapolated VPE value
        # ----------------------------------------------------------------------
        ThresholdBin = Stats.computePercentile(Cdf, 60)
        PerfInfoSer["ExtVpe"] = 5.33 * Stats.computeOverbound(Sigmas, ThresholdBin)

        # Compute continuity risk
        # ----------------------------------------------------------------------
        PerfInfoSer["ContRisk"] = PerfInfoSer["ContEvent"] / PerfInfoSer["Avail"]

        # Compute Availability Percentage
        # ----------------------------------------------------------------------
        PerfInfoSer["Avail"] = 100 * PerfInfoSer["Avail"] / PerfInfoSer["SamSol"]

    # If there are no samples available
    elif PerfInfoSer["Avail"] == 0:
        PerfInfoSer["HpeRms"] = 0.0
        PerfInfoSer["VpeRms"] = 0.0
        PerfInfoSer["Hpe95"] = 0.0
        PerfInfoSer["Vpe95"] = 0.0
        PerfInfoSer["ExtVpe"] = 0.0
        PerfInfoSer["ContRisk"] = 0.0
        PerfInfoSer["Avail"] = 0.0


def computeVpeHist(fhist, PerfLPV200, VpeHistInfo):
    # Initialize internal variables
    BinId = 0
    HistRes = GnssConstants.HIST_RES

    # Sort LPV200 VPE histogram
    SortedHist = OrderedDict({})
    for key in sorted(PerfLPV200["VpeHist"].keys()):
        SortedHist[key] = PerfLPV200["VpeHist"][key]

    # Loop over the bins in VpeHist
    for Bin, Samples in SortedHist.items():
        # Compute VPE histogram statistics
        VpeHistInfo["BinId"] = BinId
        VpeHistInfo["BinNumSam"] = Samples
        VpeHistInfo["BinFreq"] = VpeHistInfo["BinNumSam"] / (PerfLPV200["SamSol"] - PerfLPV200["NotAvail"])
        VpeHistInfo["BinMin"] = Bin
        VpeHistInfo["BinMax"] = (Bin // HistRes + 1) * HistRes
        # Update Bin ID
        BinId = BinId + 1
        # Generate output file
        generateHistFile(fhist, VpeHistInfo)
