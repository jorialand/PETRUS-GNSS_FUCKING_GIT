#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Preprocessing.py:
# This is the Preprocessing Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           Preprocessing.py
#  Date(YY/MM/DD): 01/02/21
#
#   Author: GNSS Academy
#   Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################


# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON.Utils import *
from COMMON.Iono import computeIonoMappingFunction

from InputOutput import RcvrIdx, ObsIdx, REJECTION_CAUSE
from InputOutput import FLAG, VALUE, TH, CSNEPOCHS
import numpy as np

# Preprocessing internal functions
#-----------------------------------------------------------------------


def init_output(ObsInfo, PreproObsInfo):
    """
    Fills PreproObsInfo from ObsInfo data.

    # TODO TESTING init_output: I should get a dummy ObsInfo data, and check that PreproObsInfo output is as expected.
    :param ObsInfo:
    :param PreproObsInfo:
    :return:
    """
    # Loop over satellites
    for SatObs in ObsInfo:
        # Initialize output info
        SatPreproObsInfo = {
            "Sod": 0.0,  # Second of day
            "Doy": 0,  # Day of year
            "Elevation": 0.0,  # Elevation
            "Azimuth": 0.0,  # Azimuth
            "C1": 0.0,  # GPS L1C/A pseudorange
            "P1": 0.0,  # GPS L1P pseudorange
            "L1": 0.0,  # GPS L1 carrier phase (in cycles)
            "L1Meters": 0.0,  # GPS L1 carrier phase (in m)
            "S1": 0.0,  # GPS L1C/A C/No
            "P2": 0.0,  # GPS L2P pseudorange
            "L2": 0.0,  # GPS L2 carrier phase
            "S2": 0.0,  # GPS L2 C/No
            "SmoothC1": 0.0,  # Smoothed L1CA
            "GeomFree": 0.0,  # Geom-free in Phases
            "GeomFreePrev": 0.0,  # t-1 Geom-free in Phases
            "ValidL1": 1,  # L1 Measurement Status
            "RejectionCause": 0,  # Cause of rejection flag
            "StatusL2": 0,  # L2 Measurement Status
            "Status": 0,  # L1 Smoothing status
            "RangeRateL1": 0.0,  # L1 Code Rate
            "RangeRateStepL1": 0.0,  # L1 Code Rate Step
            "PhaseRateL1": 0.0,  # L1 Phase Rate
            "PhaseRateStepL1": 0.0,  # L1 Phase Rate Step
            "VtecRate": 0.0,  # VTEC Rate
            "iAATR": 0.0,  # Instantaneous AATR
            "Mpp": 0.0,  # Iono Mapping

        }  # End of SatPreproObsInfo

        # Get satellite label
        SatLabel = SatObs[ObsIdx["CONST"]] + "%02d" % int(SatObs[ObsIdx["PRN"]])

        # Prepare outputs
        # Get SoD
        SatPreproObsInfo["Sod"] = float(SatObs[ObsIdx["SOD"]])
        # Get DoY
        SatPreproObsInfo["Doy"] = int(SatObs[ObsIdx["DOY"]])
        # Get Elevation
        SatPreproObsInfo["Elevation"] = float(SatObs[ObsIdx["ELEV"]])
        SatPreproObsInfo["Azimuth"] = float(SatObs[ObsIdx["AZIM"]])
        SatPreproObsInfo["C1"] = float(SatObs[ObsIdx["C1"]])
        # TODO WHAT IS P1¿?
        SatPreproObsInfo["L1"] = float(SatObs[ObsIdx["L1"]])
        SatPreproObsInfo["L1Meters"] = Const.GPS_L1_WAVE * float(SatObs[ObsIdx["L1"]])
        SatPreproObsInfo["S1"] = float(SatObs[ObsIdx["S1"]])
        SatPreproObsInfo["P2"] = float(SatObs[ObsIdx["P2"]])
        SatPreproObsInfo["L2"] = float(SatObs[ObsIdx["L2"]])
        SatPreproObsInfo["S2"] = float(SatObs[ObsIdx["S2"]])

        # Prepare output for the satellite
        PreproObsInfo[SatLabel] = SatPreproObsInfo

def reject_sats_lower_elev(Conf, PreproObsInfo, n_sats_to_reject):
    """
    Invalidates the n_sats_to_reject satellites with poorer visibility (i.e. elevation angle) from PreproObsInfo.
    :param Conf:
    :param PreproObsInfo: preprocessed observations
    :param n_sats_to_reject: number of sats to be rejected. (int>0)
    :return:
    """
    # Checks
    if n_sats_to_reject <= 0:
        return

    sats_elevation = {sat: sat_info['Elevation'] for (sat, sat_info) in PreproObsInfo.items()}
    for i in range(n_sats_to_reject):
        # Get the sat
        sat_to_reject = min(sats_elevation, key=sats_elevation.get)
        # Reject it
        set_sat_valid(sat_to_reject, False, REJECTION_CAUSE['NCHANNELS_GPS'], PreproObsInfo)

        del sats_elevation[sat_to_reject]

def runPreProcMeas(Conf, Rcvr, ObsInfo, PrevPreproObsInfo):
    # Purpose: preprocess GNSS raw measurements from OBS file
    #          and generate PREPRO OBS file with the cleaned,
    #          smoothed measurements

    #          More in detail, this function handles:

    #          * Measurements cleaning and validation and exclusion due to different
    #          criteria as follows:
    #             - Minimum Masking angle
    #             - Maximum Number of channels
    #             - Minimum Carrier-To-Noise Ratio (CN0)
    #             - Pseudo-Range Output of Range
    #             - Maximum Pseudo-Range Step
    #             - Maximum Pseudo-Range Rate
    #             - Maximum Carrier Phase Increase
    #             - Maximum Carrier Phase Increase Rate
    #             - Data Gaps checks and handling
    #             - Cycle Slips detection

    #         * Filtering/Smoothing of Code-Phase Measurements with a Hatch filter

    # Parameters
    # ==========
    # Conf: dict
    #         Configuration dictionary
    # Rcvr: list
    #         Receiver information: position, masking angle...
    # ObsInfo: list
    #         OBS info for current epoch
    #         ObsInfo[1][1] is the second field of the
    #         second satellite
    # PrevPreproObsInfo: dict
    #         Preprocessed observations for previous epoch per sat
    #         PrevPreproObsInfo["G01"]["C1"]

    # Returns
    # =======
    # PreproObsInfo: dict
    #         Preprocessed observations for current epoch per sat
    #         PreproObsInfo["G01"]["C1"]

    # Initialize output
    delta_t = 1
    gap_counter = {}
    rest_hatch_filter = {}

    PreproObsInfo = OrderedDict({})
    init_output(ObsInfo, PreproObsInfo)

    # Limit the satellites to the Number of Channels
    # [T2.0 CHANNELS][PETRUS-PPVE-REQ-010]
    # ----------------------------------------------------------
    n_sats_to_reject = len(PreproObsInfo.keys()) - Conf['NCHANNELS_GPS']
    if n_sats_to_reject > 0:
        reject_sats_lower_elev(Conf, PreproObsInfo, n_sats_to_reject)

    # Loop over all satellites
    for prn in range(1, Const.MAX_NUM_SATS_CONSTEL + 1):
        sat_label = get_sat_label(prn)

        # Checks
        if not sat_label in PreproObsInfo:
            continue
        if not PreproObsInfo[sat_label]['ValidL1']:
            continue

        # Reject satellites whose Elevation is lower than the Masking Angle.
        # [T2.1 ELEVATION][NO REQUIREMENT ASSOCIATED]
        # ------------------------------------------------------------------------------
        if PreproObsInfo[sat_label]["Elevation"] < Conf["RCVR_MASK"]:
            set_sat_valid(sat_label, False, REJECTION_CAUSE['MASKANGLE'], PreproObsInfo)

        # Measurement quality monitoring (SNR, PSR OUTRNG)
        # ------------------------------------------------------------------------------

        # Check Signal To Noise Ratio in front of Minimum by configuration (if activated)
        # [T2.2 SNR CHECK][PETRUS-PPVE-REQ-020]
        # ------------------------------------------------------------------------------
        check_SNR, SNR_threshold = Conf["MIN_CNR"][0], Conf["MIN_CNR"][1]
        if check_SNR:
            if PreproObsInfo[sat_label]["S1"] < SNR_threshold:
                set_sat_valid(sat_label, False, REJECTION_CAUSE['MIN_CNR'], PreproObsInfo)

        # Check Pseudo-ranges Out-of-Range in front of Maximum by configuration (if activated)
        # [T2.3 OUT-OF-RANGE][NO REQUIREMENT ASSOCIATED]
        # ------------------------------------------------------------------------------
        check_PSR, PSR_threshold = Conf["MAX_PSR_OUTRNG"][0], Conf["MAX_PSR_OUTRNG"][1]
        if check_PSR:
            if PreproObsInfo[sat_label]["C1"] > PSR_threshold:
                set_sat_valid(sat_label, False, REJECTION_CAUSE['MAX_PSR_OUTRNG'], PreproObsInfo)

        # [T2.4 GAPS][PETRUS-PPVE-REQ-090 FUN]
        if delta_t > Conf['SAMPLING_RATE']:
            gap_counter[prn] = delta_t

            if gap_counter[prn] > Conf['HATCH_GAP_TH']:
                rest_hatch_filter[prn] = True

    # Save current epoch for next iteration
    update_previous_meas(PreproObsInfo, PrevPreproObsInfo)

    return PreproObsInfo


def update_previous_meas(PreproObsInfo, PrevPreproObsInfo):
    """
    Updates PrevPreproObsInfo with current PreproObsInfo state.

    :param PreproObsInfo:
    :param PrevPreproObsInfo:
    :return:
    """
    for prn in PreproObsInfo:
        # Update carrier phase in L1
        PrevPreproObsInfo[prn]['L1_n_3'] = PrevPreproObsInfo[prn]['L1_n_2']
        PrevPreproObsInfo[prn]['L1_n_2'] = PrevPreproObsInfo[prn]['L1_n_1']
        PrevPreproObsInfo[prn]['L1_n_1'] = PreproObsInfo[prn]['L1']

        # Update epoch
        PrevPreproObsInfo[prn]['t_n_3'] = PrevPreproObsInfo[prn]['t_n_2']
        PrevPreproObsInfo[prn]['t_n_2'] = PrevPreproObsInfo[prn]['t_n_1']
        PrevPreproObsInfo[prn]['t_n_1'] = PreproObsInfo[prn]['Sod']

        # CS detection stuff...
        pass

        # Data gap detection & Smoothing
        PrevPreproObsInfo[prn]['PrevEpoch'] = PreproObsInfo[prn]['Sod']
        PrevPreproObsInfo[prn]['PrevL1'] = PreproObsInfo[prn]['L1']

        # Other stuff...
        PrevPreproObsInfo[prn]['PrevRej'] = PreproObsInfo[prn]['RejectionCause']





########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################
