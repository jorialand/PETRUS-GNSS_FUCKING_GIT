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
# ----------------------------------------------------------------------
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
# -----------------------------------------------------------------------

def check_CS_activated(Conf):
    return Conf['MIN_NCS_TH'][0] == 1

def check_PSR_activated(Conf):
    return Conf["MAX_PSR_OUTRNG"][0]

def check_SNR_activated(Conf):
    return Conf["MIN_CNR"][0] == 1

def detect_cycle_slip(meas_CS, cs_threshold, meth='TOD'):
    """
    Methods available:
        - 'TOD': Third Order Difference algorithm
    :param meas_CS: measurements used for cycle slip detection
    :param cs_threshold:
    :return:
    """
    # Checks
    if not all([meas_CS['t_n_1'], meas_CS['t_n_2'], meas_CS['t_n_3']]):
        return False

    # Previous measurements instants
    t1 = meas_CS['t_n_0'] - meas_CS['t_n_1']
    t2 = meas_CS['t_n_1'] - meas_CS['t_n_2']
    t3 = meas_CS['t_n_2'] - meas_CS['t_n_3']

    # Residuals equation factors
    # try:
    R1 = float((t1 + t2) * (t1 + t2 + t3)) / (t2 * (t2 + t3))
    R2 = float(-t1 * (t1 + t2 + t3)) / (t2 * t3)
    R3 = float(t1 * (t1 + t2)) / ((t2 + t3) * t3)
    # except ZeroDivisionError:
    #     # First iterations, may end up with t1, t2, t3 == 0
    #     return False

    # Compute TOD residuals
    CS_residual = abs(meas_CS['L1'] - R1 * meas_CS['L1_n_1'] - R2 * meas_CS['L1_n_2']  - R3 * meas_CS['L1_n_3'])

    return CS_residual > cs_threshold

def get_CS_meas(PreproObsInfo, PrevPreproObsInfo, sat_label):
    meas_CS = {}

    # L1
    meas_CS['L1'] = PreproObsInfo[sat_label]['L1']
    meas_CS['L1_n_1'] = PrevPreproObsInfo[sat_label]['L1_n_1']
    meas_CS['L1_n_2'] = PrevPreproObsInfo[sat_label]['L1_n_2']
    meas_CS['L1_n_3'] = PrevPreproObsInfo[sat_label]['L1_n_3']

    # Epoch
    meas_CS['t_n_0'] = PreproObsInfo[sat_label]['Sod']
    meas_CS['t_n_1'] = PrevPreproObsInfo[sat_label]['t_n_1']
    meas_CS['t_n_2'] = PrevPreproObsInfo[sat_label]['t_n_2']
    meas_CS['t_n_3'] = PrevPreproObsInfo[sat_label]['t_n_3']

    return meas_CS

def get_CS_threshold(Conf):
    return float(Conf['MIN_NCS_TH'][2])

def get_CS_tolerance(Conf):
    return float(Conf['MIN_NCS_TH'][1])

def get_PSR_threshold(Conf):
    return float(Conf["MAX_PSR_OUTRNG"][1])

def get_SNR_threshold(Conf):
    return float(Conf["MIN_CNR"][1])

def init_output(ObsInfo, PreproObsInfo):
    """
    Fills PreproObsInfo from ObsInfo data.

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
        # WHAT IS P1¿? TODO
        SatPreproObsInfo["L1"] = float(SatObs[ObsIdx["L1"]])
        SatPreproObsInfo["L1Meters"] = Const.GPS_L1_WAVE * float(SatObs[ObsIdx["L1"]])
        SatPreproObsInfo["S1"] = float(SatObs[ObsIdx["S1"]])
        SatPreproObsInfo["P2"] = float(SatObs[ObsIdx["P2"]])
        SatPreproObsInfo["L2"] = float(SatObs[ObsIdx["L2"]])
        SatPreproObsInfo["S2"] = float(SatObs[ObsIdx["S2"]])

        # Prepare output for the satellite
        PreproObsInfo[SatLabel] = SatPreproObsInfo

def max_consecutive_CS(PrevPreproObsInfo, sat_label, Conf):
    return sum(PrevPreproObsInfo[sat_label]['CsBuff']) == get_CS_threshold(Conf)

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

    # Build a dict with only PRN vs Elevation
    sats_elevation = {sat: sat_info['Elevation'] for (sat, sat_info) in PreproObsInfo.items()}

    # Reject satellites
    for i in range(n_sats_to_reject):
        # Get the sat
        sat_to_reject = min(sats_elevation, key=sats_elevation.get)
        # Reject it
        set_sat_valid(sat_to_reject, False, REJECTION_CAUSE['NCHANNELS_GPS'], PreproObsInfo)

        del sats_elevation[sat_to_reject]

def reset_CS_meas(PreproObsInfo, PrevPreproObsInfo, sat_label):
    # L1
    PrevPreproObsInfo[sat_label]['L1_n_1'] = PreproObsInfo[sat_label]['L1']
    PrevPreproObsInfo[sat_label]['L1_n_2'] = 0.
    PrevPreproObsInfo[sat_label]['L1_n_3'] = 0.

    # Epoch
    PrevPreproObsInfo[sat_label]['t_n_1'] = PreproObsInfo[sat_label]['Sod']
    PrevPreproObsInfo[sat_label]['t_n_2'] = 0.
    PrevPreproObsInfo[sat_label]['t_n_3'] = 0.

def runPreProcMeas(Conf, Rcvr, ObsInfo, PrevPreproObsInfo):
    # Purpose: preprocess GNSS raw measurements from OBS file
    #          and generate PREPRO OBS file with the cleaned,
    #          smoothed measurements

    #          More in detail, this function handles:

    #          * Measurements cleaning and validation and exclusion due to different
    #          criteria as follows:
    #             -- Minimum Masking angle
    #             -- Maximum Number of channels
    #             -- Minimum Carrier-To-Noise Ratio (CN0)
    #             -- Pseudo-Range Output of Range
    #             - Maximum Pseudo-Range Step
    #             - Maximum Pseudo-Range Rate
    #             - Maximum Carrier Phase Increase
    #             - Maximum Carrier Phase Increase Rate
    #             -- Data Gaps checks and handling
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
    PreproObsInfo = OrderedDict({})
    init_output(ObsInfo, PreproObsInfo)

    epoch = int(ObsInfo[0][0])
    # Print epoch (if activated)
    if TESTING and False:
        if epoch % 3600 == 0:
            print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch))
    # Limit simulation time (if activated)
    limit_time = 7200
    if TESTING and False:
        if epoch > limit_time:
            print('[TESTING][runPreProcMeas]' + ' SIMULATION TIME LIMIT REACHED (' + str(limit_time) + ')')
            exit('[TESTING][runPreProcMeas]' + ' SIMULATION TIME LIMIT REACHED (' + str(limit_time) + ')')

    # Globals
    reset_hatch_filter_vs_prn = \
        {get_sat_label(int(Sat[ObsIdx["PRN"]])): False for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}  # Container
    # to store the delta_t for all satellites

    # Limit the satellites to the Number of Channels
    # [T2.0 CHANNELS][PETRUS-PPVE-REQ-010]
    # ----------------------------------------------------------
    n_sats_to_reject = len(PreproObsInfo.keys()) - Conf['NCHANNELS_GPS']
    if n_sats_to_reject > 0:
        reject_sats_lower_elev(Conf, PreproObsInfo, n_sats_to_reject)

    # Loop over all satellites
    for prn in range(1, Const.MAX_NUM_SATS_CONSTEL + 1):
        sat_label = get_sat_label(prn)
        # reset_hatch_filter_vs_prn[sat_label] = 0

        # Checks
        if not sat_label in PreproObsInfo:
            continue
        if not PreproObsInfo[sat_label]['ValidL1']:
            continue
        # ---- From here, only sats within max channels number.
        # ------------------------------------------------------------------------------

        # Reject satellites whose Elevation is lower than the Masking Angle
        # [T2.1 ELEVATION][NO REQUIREMENT ASSOCIATED]
        # ------------------------------------------------------------------------------
        if PreproObsInfo[sat_label]["Elevation"] < Rcvr[RcvrIdx['MASK']]:
            set_sat_valid(sat_label, False, REJECTION_CAUSE['MASKANGLE'], PreproObsInfo)
            continue
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ------------------------------------------------------------------------------

        # Measurement quality monitoring (SNR, PSR OUTRNG)
        # ------------------------------------------------------------------------------

        # Check Signal To Noise Ratio in front of Minimum by configuration (if activated)
        # [T2.2 SNR CHECK][PETRUS-PPVE-REQ-020]
        # ------------------------------------------------------------------------------
        if check_SNR_activated(Conf):
            if PreproObsInfo[sat_label]["S1"] < get_SNR_threshold(Conf):
                set_sat_valid(sat_label, False, REJECTION_CAUSE['MIN_CNR'], PreproObsInfo)
                continue
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ------------------------------------------------------------------------------

        # Check Pseudo-ranges Out-of-Range in front of Maximum by configuration (if activated)
        # [T2.3 OUT-OF-RANGE][NO REQUIREMENT ASSOCIATED]
        # ------------------------------------------------------------------------------

        if check_PSR_activated(Conf):
            if PreproObsInfo[sat_label]["C1"] > get_PSR_threshold(Conf):
                set_sat_valid(sat_label, False, REJECTION_CAUSE['MAX_PSR_OUTRNG'], PreproObsInfo)
                continue
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- & max PSR
        # ------------------------------------------------------------------------------

        # # Check GAPS
        # # [T2.4 GAPS][PETRUS-PPVE-REQ-090]
        # NOTICE The delta_t should only be calculated with satellites above mask angle (aka visible sats)
        # NOTICE The non visibility period gaps [sat visible] - [sat not visible] - (gap) - [sat visible] are not
        #        considered as data gaps

        if not PrevPreproObsInfo[sat_label]['PrevEpoch']:
            # First satellite occurrence
            delta_t = Conf['SAMPLING_RATE']
        else:
            delta_t = PreproObsInfo[sat_label]['Sod'] - PrevPreproObsInfo[sat_label]['PrevEpoch']

        if delta_t > Conf['SAMPLING_RATE']:
            if delta_t > Conf['HATCH_GAP_TH']:
                reset_hatch_filter_vs_prn[sat_label] = True

                # Non-visibility periods are not gaps
                if PrevPreproObsInfo[sat_label]['PrevRej'] != REJECTION_CAUSE['MASKANGLE']:
                    set_sat_valid(sat_label, False, REJECTION_CAUSE['DATA_GAP'], PreproObsInfo)

                    if TESTING and False:
                        print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + ' Satellite ' + sat_label + \
                              ' Hatch filter reset (gap=' + "%.2f" % delta_t + ')')
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- & max PSR
        # ---- & GAPS clean
        # ------------------------------------------------------------------------------

        # Check CYCLE SLIPS (if activated)
        # [T2.5 CYCLE SLIPS][PETRUS-PPVE-REQ-080]
        # ------------------------------------------------------------------------------

        if check_CS_activated(Conf):
            if not reset_hatch_filter_vs_prn[sat_label]:
                CS_flag = detect_cycle_slip(get_CS_meas(PreproObsInfo, PrevPreproObsInfo, sat_label),
                                            get_CS_tolerance(Conf), meth='TOD')
                if CS_flag:
                    set_sat_valid(sat_label, False, REJECTION_CAUSE['CYCLE_SLIP'], PreproObsInfo)

                    if TESTING and False:
                        print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + \
                              ' Satellite ' + sat_label + '(CS)')

                CSBuff_update(CS_flag, PrevPreproObsInfo, sat_label)

                # Reset Hatch filter
                if max_consecutive_CS(PrevPreproObsInfo, sat_label, Conf):
                    reset_hatch_filter_vs_prn[sat_label] = True

                    if TESTING and True:
                        print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + \
                              ' Satellite ' + sat_label + ' Hatch filter reset (CS)')
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- & max PSR
        # ---- & GAPS clean
        # ---- & CS clean
        # ------------------------------------------------------------------------------

    # Save VALID current epoch meas for next epoch
    Ksmooth = None
    update_previous_meas(PreproObsInfo, PrevPreproObsInfo,
                         reset_hatch_filter_vs_prn, Conf, Ksmooth)

    return PreproObsInfo

def CSBuff_reset(PrevPreproObsInfo, sat_label, Conf):
    PrevPreproObsInfo[sat_label]["CsBuff"] = [0] * int(Conf["MIN_NCS_TH"][CSNEPOCHS])

def CSBuff_update(CS_flag, PrevPreproObsInfo, sat_label):
    if CS_flag:
        CSBuff_notify_CS_event(PrevPreproObsInfo, sat_label)
    else:
        CSBuff_notify_no_CS_event(PrevPreproObsInfo, sat_label)
    CSBuff_move_idx(PrevPreproObsInfo, sat_label)

def CSBuff_move_idx(PrevPreproObsInfo, sat_label):
    PrevPreproObsInfo[sat_label]["CsIdx"] += 1
    PrevPreproObsInfo[sat_label]["CsIdx"] %= len(PrevPreproObsInfo[sat_label]["CsBuff"])

def CSBuff_notify_no_CS_event(PrevPreproObsInfo, sat_label):
    PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = 0

def CSBuff_notify_CS_event(PrevPreproObsInfo, sat_label):
    PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = 1

def update_CS_buffer(PrevPreproObsInfo, sat_label, CS_flag):
    # Update CS Buffer, depending on Cycle Slip flag
    if CS_flag:
        PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = 1
    else:
        PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = 0

    # Increment CS Buffer current index circular counter
    PrevPreproObsInfo[sat_label]["CsIdx"] += 1
    PrevPreproObsInfo[sat_label]["CsIdx"] %= len(PrevPreproObsInfo[sat_label]["CsBuff"])

def update_previous_meas(PreproObsInfo, PrevPreproObsInfo, reset_hatch_filter_vs_prn, Conf, Ksmooth):
    """
    Updates PrevPreproObsInfo with current PreproObsInfo state (only VALID sats).

    :param PreproObsInfo:
    :param PrevPreproObsInfo:
    :return:
    """
    for sat_label in PreproObsInfo:
        PrevPreproObsInfo[sat_label]['PrevRej'] = PreproObsInfo[sat_label]['RejectionCause']

        # Data gap detection & Smoothing
        # Only valid measurements or when resetting the Hatch filter account for previous valid epoch
        if PreproObsInfo[sat_label]['ValidL1'] or reset_hatch_filter_vs_prn[sat_label]:
            PrevPreproObsInfo[sat_label]['PrevEpoch'] = PreproObsInfo[sat_label]['Sod']

        # CS detection
        if PreproObsInfo[sat_label]['ValidL1']:
            # Update carrier phase in L1
            PrevPreproObsInfo[sat_label]['L1_n_3'] = PrevPreproObsInfo[sat_label]['L1_n_2']
            PrevPreproObsInfo[sat_label]['L1_n_2'] = PrevPreproObsInfo[sat_label]['L1_n_1']
            PrevPreproObsInfo[sat_label]['L1_n_1'] = PreproObsInfo[sat_label]['L1']
            # Update epoch
            PrevPreproObsInfo[sat_label]['t_n_3'] = PrevPreproObsInfo[sat_label]['t_n_2']
            PrevPreproObsInfo[sat_label]['t_n_2'] = PrevPreproObsInfo[sat_label]['t_n_1']
            PrevPreproObsInfo[sat_label]['t_n_1'] = PreproObsInfo[sat_label]['Sod']

        # TAL
        if reset_hatch_filter_vs_prn[sat_label]:
            reset_CS_meas(PreproObsInfo, PrevPreproObsInfo, sat_label)
            CSBuff_reset(PrevPreproObsInfo, sat_label, Conf)

        # Valid measurements haven't reset the Hatch filter
        if PreproObsInfo[sat_label]['ValidL1']:
            PrevPreproObsInfo[sat_label]['ResetHatchFilter'] = 0

        # Hatch filter
        PrevPreproObsInfo[sat_label]["PrevL1"] = PreproObsInfo[sat_label]["L1Meters"]
        PrevPreproObsInfo[sat_label]["Ksmooth"] = Ksmooth
        PrevPreproObsInfo[sat_label]["PrevSmoothC1"] = PreproObsInfo[sat_label]["SmoothC1"]

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################
