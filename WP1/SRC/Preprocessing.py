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
import os
import sys

# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON.Utils import *
# from COMMON.Iono import computeIonoMappingFunction

from InputOutput import RcvrIdx, ObsIdx, REJECTION_CAUSE
from InputOutput import CSNEPOCHS


# Preprocessing internal functions
# -----------------------------------------------------------------------
def check_CS_activated(Conf):
    return Conf['MIN_NCS_TH'][0] == 1

def check_phase_rate_activated(Conf):
    return Conf['MAX_PHASE_RATE'][0]

def check_PSR_activated(Conf):
    return Conf["MAX_PSR_OUTRNG"][0]

def check_SNR_activated(Conf):
    return Conf["MIN_CNR"][0] == 1

def compute_phase_rate(PreproObsInfo, PrevPreproObsInfo, sat_label, delta_t):
    return Const.GPS_L1_WAVE * (PreproObsInfo[sat_label]['L1'] - PrevPreproObsInfo[sat_label]['PrevL1']) / delta_t

def detect_cycle_slip(meas_CS, cs_threshold, meth='TOD'):
    """
    Methods available:
        - 'TOD': Third Order Difference algorithm
    :param meth: method used for detecting the cycle slip event
    :param meas_CS: measurements used for cycle slip detection
    :param cs_threshold:
    :return:
    """
    if meth == 'TOD':
        # Checks
        # if not all([meas_CS['t_n_1'], meas_CS['t_n_2'], meas_CS['t_n_3']]):
        # if not meas_CS['t_n_3']:
        #     return False

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
        CS_residual = abs(meas_CS['L1'] - R1 * meas_CS['L1_n_1'] - R2 * meas_CS['L1_n_2'] - R3 * meas_CS['L1_n_3'])

        return CS_residual > cs_threshold
    else:
        exit('[runPreProcMeas][detect_cycle_slip]' + 'Not available method for cycle slip detection (' + meth + ')')

def get_CS_meas(PreproObs, PrevPreproObsInfo, sat_label):
    meas_CS = {'L1': PreproObs['L1'],
               'L1_n_1': PrevPreproObsInfo[sat_label]['L1_n_1'],
               'L1_n_2': PrevPreproObsInfo[sat_label]['L1_n_2'],
               'L1_n_3': PrevPreproObsInfo[sat_label]['L1_n_3'],
               't_n_0': PreproObs['Sod'],
               't_n_1': PrevPreproObsInfo[sat_label]['t_n_1'],
               't_n_2': PrevPreproObsInfo[sat_label]['t_n_2'],
               't_n_3': PrevPreproObsInfo[sat_label]['t_n_3']}
    return meas_CS

def get_CS_threshold(Conf):
    return float(Conf['MIN_NCS_TH'][2])

def get_CS_tolerance(Conf):
    return float(Conf['MIN_NCS_TH'][1])

def get_elev_cut_value_nchannels(PreproObsInfo, Conf):
    """
    Get the elevation value starting from which to reject satellites due to maximum channels limit.
    :param PreproObsInfo:
    :param Conf:
    :return: elevation
    """
    elevation_cut_value_nchannels = 0.
    n_sats_to_reject = len(PreproObsInfo) - int(Conf['NCHANNELS_GPS'])

    if n_sats_to_reject > 0:
        elevations_list = []
        for sat_label, PreproObs in PreproObsInfo.items():
            elevations_list.append(PreproObs['Elevation'])
        elevations_list = sorted(elevations_list)

        elevation_cut_value_nchannels = elevations_list[n_sats_to_reject]
    return elevation_cut_value_nchannels

def get_hatch_gap_threshold(Conf):
    return Conf['HATCH_GAP_TH']

def get_hatch_time_threshold(Conf):
    return Conf['HATCH_TIME']

def get_phase_rate_threshold(Conf):
    return Conf['MAX_PHASE_RATE'][1]

def get_pred_C1(PreproObsInfo, PrevPreproObsInfo, sat_label):
    C1 = PreproObsInfo[sat_label]['C1']
    pred_C1_prev = PrevPreproObsInfo[sat_label]['PrevSmoothC1']
    L1 = PreproObsInfo[sat_label]['L1']
    L1_prev = PrevPreproObsInfo[sat_label]['PrevL1']
    return C1, pred_C1_prev, L1, L1_prev

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
        SatPreproObsInfo["Sod"] = float(SatObs[ObsIdx["SOD"]])
        SatPreproObsInfo["Doy"] = int(SatObs[ObsIdx["DOY"]])
        SatPreproObsInfo["Elevation"] = float(SatObs[ObsIdx["ELEV"]])
        SatPreproObsInfo["Azimuth"] = float(SatObs[ObsIdx["AZIM"]])
        SatPreproObsInfo["C1"] = float(SatObs[ObsIdx["C1"]])
        SatPreproObsInfo["L1"] = float(SatObs[ObsIdx["L1"]])
        SatPreproObsInfo["L1Meters"] = Const.GPS_L1_WAVE * float(SatObs[ObsIdx["L1"]])
        SatPreproObsInfo["S1"] = float(SatObs[ObsIdx["S1"]])
        SatPreproObsInfo["L2"] = float(SatObs[ObsIdx["L2"]])
        # Not needed P2 & S2
        # SatPreproObsInfo["P2"] = float(SatObs[ObsIdx["P2"]])
        # SatPreproObsInfo["S2"] = float(SatObs[ObsIdx["S2"]])

        # Prepare output for the satellite
        PreproObsInfo[SatLabel] = SatPreproObsInfo

def max_consecutive_CS(PrevPreproObsInfo, sat_label, Conf):
    return sum(PrevPreproObsInfo[sat_label]['CsBuff']) == get_CS_threshold(Conf), sum(PrevPreproObsInfo[sat_label]['CsBuff'])

def reset_CS_meas(PreproObs, PrevPreproObsInfo, sat_label, hard_reset=False):
    # L1
    PrevPreproObsInfo[sat_label]['L1_n_1'] = PreproObs['L1']
    PrevPreproObsInfo[sat_label]['L1_n_2'] = 0.
    PrevPreproObsInfo[sat_label]['L1_n_3'] = 0.

    # Epoch
    PrevPreproObsInfo[sat_label]['t_n_1'] = PreproObs['Sod']
    PrevPreproObsInfo[sat_label]['t_n_2'] = 0.
    PrevPreproObsInfo[sat_label]['t_n_3'] = 0.

    if hard_reset:
        PrevPreproObsInfo[sat_label]['L1_n_1'] = 0.
        PrevPreproObsInfo[sat_label]['t_n_1'] = 0.

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

def CSBuff_notify_CS_event(PrevPreproObsInfo, sat_label, CS_flag):
    PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = CS_flag

def reset_hatch_filter(PreproObs, PrevPreproObsInfo, sat_label, Conf):
    epoch = PreproObs['Sod']
    # Reset gap counter
    PrevPreproObsInfo[sat_label]['gap_counter'] = 0.
    # Reset smoothing artifacts
    PrevPreproObsInfo[sat_label]['k_smooth'] = 1
    PreproObs['SmoothC1'] = PreproObs['C1']
    PrevPreproObsInfo[sat_label]['PrevSmoothC1'] = PreproObs['SmoothC1']
    # Reset phase measurements
    PrevPreproObsInfo[sat_label]['PrevL1'] = PreproObs['L1']
    PrevPreproObsInfo[sat_label]['PrevEpoch'] = epoch
    # Reset rates
    PrevPreproObsInfo[sat_label]['PrevRangeRateL1'] = -9999.99
    PrevPreproObsInfo[sat_label]['PrevPhaseRateL1'] = -9999.99
    # More smoothing artifacts
    PrevPreproObsInfo[sat_label]['reset_hatch_filter'] = False
    # Update smoothing status (0 means smoothing isn't done)
    PreproObs['Status'] = 0
    # Reset CS detection artifacts
    reset_CS_meas(PreproObs, PrevPreproObsInfo, sat_label, hard_reset=True)
    CSBuff_reset(PrevPreproObsInfo, sat_label, Conf)

def update_CS_buffer(PrevPreproObsInfo, sat_label, CS_flag):
    # Update CS Buffer, depending on Cycle Slip flag
    if CS_flag:
        PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = 1
    else:
        PrevPreproObsInfo[sat_label]["CsBuff"][PrevPreproObsInfo[sat_label]["CsIdx"]] = 0

    # Increment CS Buffer current index circular counter
    PrevPreproObsInfo[sat_label]["CsIdx"] += 1
    PrevPreproObsInfo[sat_label]["CsIdx"] %= len(PrevPreproObsInfo[sat_label]["CsBuff"])

# MOST RELEVANT FUNCTION HERE
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
    # Print % progress (if activated)
    if TESTING and True:
        step = epoch / 86400 * 100
        if step % 10 == 0:
            print('[TESTING][runPreProcMeas]' + ' Progress ' + '%.d' % step + '%')
    # Limit simulation time (if activated)
    limit_time = 7200
    if TESTING and False:
        if epoch > limit_time:
            print('[TESTING][runPreProcMeas]' + ' SIMULATION TIME LIMIT REACHED (' + str(limit_time) + ')')
            exit('[TESTING][runPreProcMeas]' + ' SIMULATION TIME LIMIT REACHED (' + str(limit_time) + ')')
    TESTING_PRINT = {'nchannels': True,     # Tested, none
                     'maskangle': False,    # Tested, many
                     'snr': False,          # Tested, some
                     'psr': True,           # Tested, none
                     'reset_hatch_filter': True, # Tested, many
                     'gap': False,           # Tested, some
                     'CS': True,            # Tested, some
                     }
    # Globals
    pass

    # Get the minimum elevation value. All sats below it, will be rejected due to maximum nchannels
    elev_cut_value_nchannels = get_elev_cut_value_nchannels(PreproObsInfo, Conf)

    # Loop over all satellites
    for sat_label, PreproObs in PreproObsInfo.items():
        # Limit the satellites to the Number of Channels
        # [T2.0 CHANNELS][PETRUS-PPVE-REQ-010]
        # ----------------------------------------------------------
        if PreproObs['Elevation'] < elev_cut_value_nchannels:
            set_sat_valid(PreproObs, False, REJECTION_CAUSE['NCHANNELS_GPS'])
            if TESTING and TESTING_PRINT['nchannels']:
                print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) +
                      ' Satellite ' + sat_label + ' Rejected (NCHANNELS_GPS)')
            continue
        # ---- From here, only sats within max channels number.
        # ---- TESTED WITH ./testing_nchannels.sh
        # ------------------------------------------------------------------------------

        # Reject satellites whose Elevation is lower than the Masking Angle
        # [T2.1 ELEVATION][NO REQUIREMENT ASSOCIATED]
        # ------------------------------------------------------------------------------
        if PreproObs["Elevation"] < Rcvr[RcvrIdx['MASK']]:
            set_sat_valid(PreproObs, False, REJECTION_CAUSE['MASKANGLE'])
            PrevPreproObsInfo[sat_label]['PrevRej'] = REJECTION_CAUSE['MASKANGLE']
            if TESTING and TESTING_PRINT['maskangle']:
                print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) +
                      ' Satellite ' + sat_label + ' Rejected (MASKANGLE)')
            continue
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- TESTED WITH ./testing_maskangle.sh
        # ------------------------------------------------------------------------------

        # Measurement quality monitoring (SNR, PSR OUTRNG)
        # ------------------------------------------------------------------------------

        # Check Signal To Noise Ratio in front of Minimum by configuration (if activated)
        # [T2.2 SNR CHECK][PETRUS-PPVE-REQ-020]
        # ------------------------------------------------------------------------------
        if check_SNR_activated(Conf):
            if PreproObs["S1"] < get_SNR_threshold(Conf):
                set_sat_valid(PreproObs, False, REJECTION_CAUSE['MIN_CNR'])
                PrevPreproObsInfo[sat_label]['PrevRej'] = REJECTION_CAUSE['MIN_CNR']
                if TESTING and TESTING_PRINT['snr']:
                    print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) +
                          ' Satellite ' + sat_label + ' Rejected (MIN_CNR)')
                continue
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- TESTED WITH ./testing_snr.sh
        # ------------------------------------------------------------------------------

        # Check Pseudo-ranges Out-of-Range in front of Maximum by configuration (if activated)
        # [T2.3 OUT-OF-RANGE][NO REQUIREMENT ASSOCIATED]
        # ------------------------------------------------------------------------------
        if check_PSR_activated(Conf):
            if PreproObs["C1"] > get_PSR_threshold(Conf):
                set_sat_valid(PreproObs, False, REJECTION_CAUSE['MAX_PSR_OUTRNG'])
                # PrevPreproObsInfo[sat_label]['PrevRej'] = REJECTION_CAUSE['MAX_PSR_OUTRNG'] # TODO PrevRej MAX_PSR_OUTRNG¿?
                if TESTING and TESTING_PRINT['psr']:
                    print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) +
                          ' Satellite ' + sat_label + ' Rejected (MAX_PSR)')
                continue
        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- & max PSR
        # ---- TESTED WITH ./testing_psr.sh
        # ------------------------------------------------------------------------------

        # Check GAPS
        # [T2.4 GAPS][PETRUS-PPVE-REQ-090]
        # ------------------------------------------------------------------------------
        # NOTICE The delta_t should only be calculated with satellites above mask angle (aka visible sats)
        # NOTICE The non visibility period gaps [sat visible] - [sat not visible] - (gap) - [sat visible] are not
        #        considered as data gaps
        # NOTICE First satellite occurrence doesn't account as data gap bcause delta_t<0
        # NOTICE Resuming: first, when a data gap is detected, reset hatch filter (this includes non-visibilty periods.
        #       Furthermore, raise the data gap detection flag (excluding sats coming from non-visible periods) (but do
        #       not invalidate the sat, nor reject it) TODO WHY NOT INVALIDATE & REJECT DATA GAPS¿?
        delta_t = epoch - PrevPreproObsInfo[sat_label]['PrevEpoch']
        if delta_t > Conf['SAMPLING_RATE']:
            PrevPreproObsInfo[sat_label]['gap_counter'] = delta_t
            if PrevPreproObsInfo[sat_label]['gap_counter'] > get_hatch_gap_threshold(Conf):
                PrevPreproObsInfo[sat_label]['reset_hatch_filter'] = True
                PrevPreproObsInfo[sat_label]['gap_counter'] = 0

                # NOTICE Non-visibility periods are not gaps. In case the previous rejection cause was bcause mask angle,
                #        it means that the satellite has become visible after a period of non-visibility
                if PrevPreproObsInfo[sat_label]['PrevRej'] != REJECTION_CAUSE['MASKANGLE']:
                    set_sat_valid(PreproObs, True, REJECTION_CAUSE['DATA_GAP'])  # TODO WHY NOT INVALIDATE & REJECT DATA GAPS¿?
                    if TESTING and TESTING_PRINT['gap']:
                        print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) + ' Satellite ' + sat_label +
                              ' Rejected* (gap=' + "%.2f" % delta_t + ')')
                else:
                    if TESTING and TESTING_PRINT['reset_hatch_filter']:
                        print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) + ' Satellite ' + sat_label +
                              ' Hatch filter reset (non-visibility gap=' + "%.2f" % delta_t + ')')
        else:
            PrevPreproObsInfo[sat_label]['gap_counter'] = 0

        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- & max PSR
        # ---- & GAPS clean
        # ---- TESTED WITH ./testing_gaps.sh
        # ------------------------------------------------------------------------------

        # Check CYCLE SLIPS (if activated)
        # [T2.5 CYCLE SLIPS][PETRUS-PPVE-REQ-080]
        # ------------------------------------------------------------------------------
        # NOTICE If no 3 previous valid meas, just update the CS detection buffer
        #        Else compute the CS detection with TOD
        # NOTICE If CS detected, reject satellite (w/o rejection cause)
        # NOTICE This will generate sats not valid with 0 rejection cause
        # NOTICE Else just update the CS detection buffer
        # NOTICE If CS detected & full CS buffer, add rejection cause CS & reset Hatch filter
        if sat_label == 'G01':
            stop = True
            if epoch == 6480:
                stop = True
        if check_CS_activated(Conf):
            if not PrevPreproObsInfo[sat_label]['reset_hatch_filter']:
                if PrevPreproObsInfo[sat_label]['t_n_3']:
                    CS_flag = detect_cycle_slip(get_CS_meas(PreproObs, PrevPreproObsInfo, sat_label),
                                                get_CS_tolerance(Conf), meth='TOD')
                    CSBuff_notify_CS_event(PrevPreproObsInfo, sat_label, CS_flag)
                    # Just invalidate the satellite because of phase cycle slip (no reject till successive CS detected)
                    if CS_flag:
                        # Invalidate the satellite
                        set_sat_valid(PreproObs, False, None)
                        if TESTING and TESTING_PRINT['CS']:
                            print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) + ' Satellite ' + sat_label +
                                  ' (CS flag)')
                        # CS Buffer full¿? -> reset Hatch filter
                        CS_buffer_full, sum = max_consecutive_CS(PrevPreproObsInfo, sat_label, Conf)
                        if CS_buffer_full:
                            set_sat_valid(PreproObs, None, REJECTION_CAUSE['CYCLE_SLIP'])
                            PrevPreproObsInfo[sat_label]['reset_hatch_filter'] = True
                            if TESTING and TESTING_PRINT['CS']:
                                print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) + ' Satellite ' + sat_label +
                                      ' (Full CS buffer)')
                        # Update CS detection buffer and skip satellite
                        else:
                            CSBuff_move_idx(PrevPreproObsInfo, sat_label)
                            continue

                # Update CS detection buffer
                CSBuff_move_idx(PrevPreproObsInfo, sat_label)
                CS_buffer_full, sum = max_consecutive_CS(PrevPreproObsInfo, sat_label, Conf)
                # Store L1 & epoch values for TOD CS detection
                if sum == 0:
                    PrevPreproObsInfo[sat_label]['L1_n_3'] = PrevPreproObsInfo[sat_label]['L1_n_2']
                    PrevPreproObsInfo[sat_label]['L1_n_2'] = PrevPreproObsInfo[sat_label]['L1_n_1']
                    PrevPreproObsInfo[sat_label]['L1_n_1'] = PreproObs['L1']
                    PrevPreproObsInfo[sat_label]['t_n_3'] = PrevPreproObsInfo[sat_label]['t_n_2']
                    PrevPreproObsInfo[sat_label]['t_n_2'] = PrevPreproObsInfo[sat_label]['t_n_1']
                    PrevPreproObsInfo[sat_label]['t_n_1'] = epoch


        # ---- From here, only sats within max channels number
        # ---- & min mask angle
        # ---- & min SNR
        # ---- & max PSR
        # ---- & GAPS clean
        # ---- & CS clean
        # ---- TESTED WITH ./testing_cs.sh
        # ------------------------------------------------------------------------------

    #     # Hatch filter (re)initialization
    #     # [T2.6 SMOOTHING][PETRUS-PPVE-REQ-100]
    #     # -------------------------------------------------------------------------------
    #     # Reset Hatch filter
        if PrevPreproObsInfo[sat_label]['reset_hatch_filter']:
            reset_hatch_filter(PreproObs, PrevPreproObsInfo, sat_label, Conf)
            if TESTING and TESTING_PRINT['reset_hatch_filter']:
                print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) + ' Satellite ' + sat_label +
                      ' Hatch filter reset')
            continue
    #
    #     #  Perform the Code Carrier SMOOTHING with a Hatch Filter
    #     # -------------------------------------------------------------------------------
    #     # Count the number of seconds from Smoothing Start
    #     PrevPreproObsInfo[sat_label]['k_smooth'] += delta_t
    #
    #     # Update the Smoothing time
    #     # smooth_time = k_smooth if k_smooth <= HATCH_TIME
    #     # smooth_time = HATCH_TIME if k_smooth > HATCH_TIME
    #     smooth_time = \
    #         (PrevPreproObsInfo[sat_label]['k_smooth'] <= get_hatch_time_threshold(Conf)) \
    #         * PrevPreproObsInfo[sat_label]['k_smooth'] \
    #                 + (PrevPreproObsInfo[sat_label]['k_smooth'] > get_hatch_time_threshold(Conf)) \
    #         * get_hatch_time_threshold(Conf)
    #
    #     # Weighting factor of the smoothing filter
    #     alpha = delta_t / smooth_time
    #
    #     if TESTING and False and sat_label == 'G10':
    #         print('[TESTING][runPreProcMeas]' + ' epoch' + str(epoch) +
    #               ' Satellite ' + sat_label + ' alpha=' + str(alpha))
    #
    #     # if sat_label == 'G10':
    #     #     print('epoch' + str(epoch) + ' alpha=' + str(alpha))
    #     # Compute predicted smooth code
    #     C1, pred_C1_prev, L1, L1_prev = get_pred_C1(PreproObsInfo, PrevPreproObsInfo, sat_label)
    #     pred_C1 = C1 + (pred_C1_prev + (L1 - L1_prev) * Const.GPS_L1_WAVE)
    #
    #     # Compute the new Smoothed Code
    #     PreproObsInfo[sat_label]['SmoothC1'] = alpha * PreproObsInfo[sat_label]['C1'] + (1 - alpha) * pred_C1
    #
    #     # Check Phase Rate (if activated)
    #     # PreproObsInfo[sat_label]['PhaseRateL1'] = \
    #     #     compute_phase_rate(PreproObsInfo, PrevPreproObsInfo, sat_label, delta_t)
    #     #
    #     # if check_phase_rate_activated(Conf):
    #     #     if PreproObsInfo[sat_label]['PhaseRateL1'] > get_phase_rate_threshold(Conf):
    #     #         PrevPreproObsInfo[sat_label]['reset_hatch_filter'] = True
    #     #         set_sat_valid(sat_label, False, REJECTION_CAUSE['MAX_PHASE_RATE'], PreproObsInfo)
    #     #         continue
    #     # Check Phase Rate Step (if activated)
    #     pass
    #
    #     # Check Code Rate detector (if activated)
    #     pass
    #
    #     # Check Code Rate Step detector
    #     pass
    #
    #     # Update Measurement Smoothing Status and Hatch filter Convergence
    #     pass
        # Save VALID current epoch meas for next epoch
        update_previous_meas(PreproObs, PrevPreproObsInfo, sat_label)


    return PreproObsInfo

# 2nd MOST RELEVANT FUNCTION HERE
def update_previous_meas(PreproObs, PrevPreproObsInfo, sat_label):
    epoch = PreproObs['Sod']
    # PrevPreproObsInfo[sat_label]['PrevSmoothC1'] = PreproObs['SmoothC1']
    PrevPreproObsInfo[sat_label]['PrevL1'] = PreproObs['L1']
    PrevPreproObsInfo[sat_label]['PrevEpoch'] = epoch
    # PrevPreproObsInfo[sat_label]['PrevRangeRateL1'] = PreproObs['RangeRateL1']
    # PrevPreproObsInfo[sat_label]['PrevPhaseRateL1'] = PreproObs['PhaseRateL1']
    PrevPreproObsInfo[sat_label]['PrevRej'] = PreproObs['RejectionCause']

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################
