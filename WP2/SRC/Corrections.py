#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Corrections.py:
# This is the Corrections Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           Corrections.py
#  Date(YY/MM/DD): 16/02/21
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
from InputOutput import RcvrIdx, SatIdx, LosIdx
import numpy as np

# Compute the clock relativistic effect
def computeDtr(SatInfo, SatLabel):
    x = float(SatInfo[SatLabel][SatIdx['SAT-X']])
    y = float(SatInfo[SatLabel][SatIdx['SAT-Y']])
    z = float(SatInfo[SatLabel][SatIdx['SAT-Z']])
    vx = float(SatInfo[SatLabel][SatIdx['VEL-X']])
    vy = float(SatInfo[SatLabel][SatIdx['VEL-Y']])
    vz = float(SatInfo[SatLabel][SatIdx['VEL-Z']])
    return (-2 * np.dot([x,y,z], [vx,vy,vz]) ) / Const.SPEED_OF_LIGHT

# Compute the obliquity factor Fpp for UISD & UIRE
def computeObliquityFactorFpp(LosInfo, SatLabel):
    return ( 1-(
            (Const.EARTH_RADIUS * np.cos(np.deg2rad(LosInfo[SatLabel][LosIdx['ELEV']]))) /
            (Const.EARTH_RADIUS + Const.IONO_HEIGHT) )**2 )**(-0.5)

# Compute the SigmaFLT projected into the User direction
def computeSigmaFlt(SatInfo, SatLabel):
    """
    Applicable standards: MOPS-DO-229D Section A.4.5.1
    """
    sigmaFLT: float

    # Get the data from SatInfo
    RSS_flag = bool(SatInfo[SatLabel][SatIdx['RSS']])
    sigma_udre = float(SatInfo[SatLabel][SatIdx['SIGMAUDRE']])
    delta_udre = float(SatInfo[SatLabel][SatIdx['DELTAUDRE']])
    epsilon_fc = float(SatInfo[SatLabel][SatIdx['EPS-FC']])
    epsilon_rrc = float(SatInfo[SatLabel][SatIdx['EPS-RRC']])
    epsilon_ltc = float(SatInfo[SatLabel][SatIdx['EPS-LTC']])
    epsilon_er = float(SatInfo[SatLabel][SatIdx['EPS-ER']])

    # Perform the SigmaFLT correction
    # Computation depending on RRS flag (Root Sum Square) from MessageType10
    if RSS_flag:
        sigmaFLT = sum(i*i for i in [sigma_udre*delta_udre, epsilon_fc, epsilon_rrc, epsilon_ltc, epsilon_er])
    else:
        sigmaFLT = ((sigma_udre * delta_udre) + epsilon_fc + epsilon_rrc + epsilon_ltc + epsilon_er)**2

    return sigmaFLT

# Compute the ionospheric corrections and error bounds UISD & UIRE
def computeUISDandUIRE(LosInfo, SatLabel):
    """
    Applicable standards: MOPS-DO-229D Section A.4.4.10.3
    """
    # Checks
    if not int(LosInfo[SatLabel][LosIdx['FLAG']]):
        sys.stderr.write("[WARNING][Corrections][computeUISDandUIRE] Not valid computation for PA.")

    # Square interpolation
    if LosInfo[SatLabel][LosIdx['INTERP']] == 0:
        # Compute GIVD @ IPP

        # Coordinates computation
        phi_pp = LosInfo[SatLabel][LosIdx['IPPLAT']] # IPP's latitude
        lambda_pp = LosInfo[SatLabel][LosIdx['IPPLON']] # IPP's longitude

        phi_1 = LosInfo[SatLabel][LosIdx['IGP_SW_LAT']] # According to FIGURE A-19
        phi_2 = LosInfo[SatLabel][LosIdx['IGP_NW_LAT']]
        lambda_1 = LosInfo[SatLabel][LosIdx['IGP_SW_LON']]
        lambda_2 = LosInfo[SatLabel][LosIdx['IGP_SE_LON']]

        delta_phi_pp = phi_pp - phi_1
        delta_lambda_pp = lambda_pp - lambda_1

        # For mid latitudes
        if phi_pp > -85. and phi_pp < 85.:
            x_pp = delta_lambda_pp / (lambda_2 - lambda_1)
            y_pp = delta_phi_pp / (phi_2 - phi_1)
        # For extreme latitudes
        else:
            # Did not found how to compute it, by comparing MOPS guidelines with PETRUS LOS file
            # Assumed that we are not treating extreme latitudes
            x_pp = delta_lambda_pp / (lambda_2 - lambda_1)
            y_pp = delta_phi_pp / (phi_2 - phi_1)
            sys.stderr.write("[WARNING][Corrections][computeUISDandUIRE] Found extreme IPP latitude.")

        # Weighting functions
        w_1 = x_pp * y_pp
        w_2 = (1-x_pp) * y_pp
        w_3 = (1-x_pp) * (1-y_pp)
        w_4 = x_pp * (1-y_pp)

        # Compute interpolated (Square) GIVD @ IPP

        GIVD = sum(W_i*GIVD_i for W_i,GIVD_i in {w_1:LosInfo[SatLabel][LosIdx['GIVD_NE']],
                                                 w_2:LosInfo[SatLabel][LosIdx['GIVD_NW']],
                                                 w_3:LosInfo[SatLabel][LosIdx['GIVD_SW']],
                                                 w_4:LosInfo[SatLabel][LosIdx['GIVD_SE']]}.items())

        # UISD
        Fpp = computeObliquityFactorFpp(LosInfo, SatLabel)
        UISD = -Fpp * GIVD

        # GIVE @ IPP (hay que interpolar papá)
        GIVE = 1.
        pass
        UIRE = np.sqrt(Fpp ** 2 * GIVE ** 2)

    # Triangle interpolation
    else:
        # interp_type represents the vertex of the square not used in the triangle interpolation
        pass






    return UISD, UIRE

# Compute thee Satellite Corrected Position and Clock applying the SBAS FLT Corrections
def correctSatPosClk(SatInfo, SatLabel):
    # Satellite position at TT corrected by SBAS LTC (Also corrected from Sagnac Effect. WGS84 ref)
    sbas_x = float(SatInfo[SatLabel][SatIdx['SAT-X']]) + float(SatInfo[SatLabel][SatIdx['LTC-X']])
    sbas_y = float(SatInfo[SatLabel][SatIdx['SAT-Y']]) + float(SatInfo[SatLabel][SatIdx['LTC-Y']])
    sbas_z = float(SatInfo[SatLabel][SatIdx['SAT-Z']]) + float(SatInfo[SatLabel][SatIdx['LTC-Z']])

    # Satellite clock corrected with SBAS FLT
    b_sbas = float(SatInfo[SatLabel][SatIdx['SAT-CLK']])  + computeDtr(SatInfo, SatLabel) - \
             float(SatInfo[SatLabel][SatIdx['TGD']]) + float(SatInfo[SatLabel][SatIdx['FC']]) + \
             float(SatInfo[SatLabel][SatIdx['LTC-B']])

    # Checks
    if 0. in [sbas_x, sbas_y, sbas_z, b_sbas]:
        sys.stderr.write("[WARNING][Corrections][correctSatPosClk] Null data when it should't")

    return (sbas_x, sbas_y, sbas_z), b_sbas


# The most important function here
def runCorrectMeas(Conf, Rcvr, PreproObsInfo, SatInfo, LosInfo):

    # Purpose: correct GNSS preprocessed measurements and compute
    #          pseudo range residuals

    #          More in detail, this function handles the following:
    #          tasks:

    #             *  Correct the satellite navigation position and clock using EGNOS Fast-Long-Term (FLT) corrections: FC and LTC.
    #             *  Estimate the Slant Ionospheric Delay (UISD) using MOPS guidelines interpolation criteria for IGP Selection
    #             *  Estimate the Slant Troposphere delay (STD) using MOPS model (ZTD) and its mapping function. 
    #             *  Correct the Pre-processed measurements from Geometrical Range, Satellite clock, ionosphere and troposphere. 
    #             *  Build the Corrected Measurements and Measurement Residuals
    #             *  Estimate all Range level Sigmas to build the Sigma UERE:
    #                   -  Estimate the SigmaUIRE from MT26 information
    #                   -  Estimate the SigmaFLT from UDRE and MT28 
    #                   -  Estimate the SigmaTRO budget in line with MOPS.
    #                   -  Estimate the SigmaAirborne budget in line with MOPS 

    #             *  Estimate the Sigma UERE budget in line with MOPS


    # Parameters
    # ==========
    # Conf: dict
    #         Configuration dictionary
    # Rcvr: list
    #         Receiver information: position, masking angle...
    # PreproObsInfo: dict
    #         Preprocessed observations for current epoch per sat
    #         PreproObsInfo["G01"]["C1"]
    # SatInfo: dict
    #         dictionary containing the split lines of the SAT file
    #         SatInfo["G01"][1] is the second field of the line
    #         containing G01 info
    # LosInfo: dict
    #         dictionary containing the split lines of the LOS file
    #         SatInfo["G01"][1] is the second field of the line
    #         containing G01 info

    # Returns
    # =======
    # CorrInfo: dict
    #         Corrected measurements for current epoch per sat
    #         CorrInfo["G01"]["CorrectedPsr"]

    # Initialize output
    CorrInfo = OrderedDict({})

    # Initialize some values
    ResSum = 0.0
    ResN = 0
    EntGpsSum = 0.0
    EntGpsN = 0

    # Loop over satellites
    for SatLabel, SatPrepro in PreproObsInfo.items():
        epoch = SatPrepro['Sod']

        # If satellite is in convergence
        if(SatPrepro["Status"] == 1):
            # Initialize output info
            SatCorrInfo = {
                "Sod": 0.0,             # Second of day
                "Doy": 0,               # Day of year
                "Elevation": 0.0,       # Elevation
                "Azimuth": 0.0,         # Azimuth
                "IppLon": 0.0,          # IPP Longitude
                "IppLat": 0.0,          # IPP Latitude
                "Flag": 1,              # 0: Not Used 1: Used for PA 2: Used for NPA
                "SatX": 0.0,            # X-Component of the Satellite Position corrected with SBAS LTC
                "SatY": 0.0,            # Y-Component of the Satellite Position corrected with SBAS LTC
                "SatZ": 0.0,            # Z-Component of the Satellite Position corrected with SBAS LTC
                "SatClk": 0.0,          # Satellite Clock corrected with SBAS FLT
                "Uisd": 0.0,            # User Ionospheric Slant Delay
                "Std": 0.0,             # Slant Tropospheric Delay
                "CorrPsr": 0.0,         # Pseudo Range corrected from delays
                "GeomRange": 0.0,       # Geometrical Range (distance between Satellite Position and Receiver Reference Position)
                "PsrResidual": 0.0,     # Pseudo Range Residual
                "RcvrClk": 0.0,         # Receiver Clock estimation
                "SigmaFlt": 0,          # Sigma of the residual error associated to the fast and long-term correction (FLT)
                "SigmaUire": 0,         # User Ionospheric Range Error Sigma
                "SigmaTropo": 0,        # Sigma of the Tropospheric error
                "SigmaAirborne": 0.0,   # Sigma Airborne Error
                "SigmaNoiseDiv": 0.0,   # Sigma of the receiver noise + divergence
                "SigmaMultipath": 0.0,  # Sigma of the receiver multipath
                "SigmaUere": 0.0,       # Sigma User Equivalent Range Error (Sigma of the total residual error associated to the satellite)
                "EntGps": 0.0,          # ENT to GPS Offset

            } # End of SatCorrInfo

            # Prepare outputs
            # Get SoD
            SatCorrInfo["Sod"] = SatPrepro["Sod"]
            # Get DoY
            SatCorrInfo["Doy"] = SatPrepro["Doy"]
            # Get Elevation
            SatCorrInfo["Elevation"] = SatPrepro["Elevation"]
            # Get Azimuth
            SatCorrInfo["Azimuth"] = SatPrepro["Azimuth"]

            # If SBAS information is available for current satellite
            # (Cannot continue if for current epoch, the current satellite is
            # not present in the Sat/Los files aka MONITORED)
            if not ((SatLabel in SatInfo) and (SatLabel in LosInfo)):
                continue

            # Get IPP Longitude
            SatCorrInfo["IppLon"] = float(LosInfo[SatLabel][LosIdx["IPPLON"]])
            # Get IPP Latitude
            SatCorrInfo["IppLat"] = float(LosInfo[SatLabel][LosIdx["IPPLAT"]])

            # NOTE: According to this, if one satellite has not IppLon or IppLat fields filled, this means it
            #       is not being monitored.

            # Check monitoring status (UDREi<12 for Service level “PA”) and (FLAG==1 for SL "PA")
            is_monitored_for_PA = (int(SatInfo[SatLabel][SatIdx['UDREI']]) < 12 ) and \
                                  (int(LosInfo[SatLabel][LosIdx['FLAG']]))
            if is_monitored_for_PA:
                # Apply the SBAS corrections to the satellite position and clock
                # [T2.1.1 SAT CORRECTION AND SIGMA FLT][PETRUS-CORR-REQ-010]
                # ----------------------------------------------------------------------
                SatPos, SatClk = correctSatPosClk(SatInfo,SatLabel)
                SatCorrInfo['SatX'] = SatPos[0]
                SatCorrInfo['SatY'] = SatPos[1]
                SatCorrInfo['SatZ'] = SatPos[2]
                SatCorrInfo['SatClk'] = SatClk

                # Compute the SigmaFLT projected into the User direction
                # [T2.1.2 SAT CORRECTION AND SIGMA FLT][PETRUS-CORR-REQ-030]
                # ----------------------------------------------------------------------
                SatCorrInfo['SigmaFlt'] = computeSigmaFlt(SatInfo, SatLabel)

                # Compute UISD and UIRE on the IPP using MOPS interpolation
                # [T2.2 UISD & UIRE][PETRUS-CORR-REQ-050][PETRUS-CORR-REQ-070]
                # ----------------------------------------------------------------------
                UISD, UIRE= computeUISDandUIRE(LosInfo, SatLabel)








                # -----------------------------------------------------------------------
                # Compute Tropospheric Mapping Function
                # TropoMpp = computeTropoMpp(Elevation)

                # SCHEMATIC

                # Compute the Slant Tropospheric Delay Error Sigma
                # SigmaTROPO = computeSigmaTROPO(TropoMpp)

                # SCHEMATIC

                # -----------------------------------------------------------------------
                # [CHALLENGING OPTION] Compute the Slant Tropospheric Delay
                # Note that STD can be read from LOS file, or
                # Challenging option using MOPS Tropo Model in Appendix A.
                # STD = computeSlantTropoDelay(RCVR[iRec].llh, Doy, TropoMpp)

                # SCHEMATIC
                # -----------------------------------------------------------------------

                # Compute User Airborne Sigma. Ref: MOPS-DO-229D Section J.2.4
                # -----------------------------------------------------------------------
                # Consider Maximum Signal Level when satellite elevation is greater
                # than Conf.ELEV_NOISE_TH=20, and Minimum Signal Level otherwise
                # SigmaAIR = computeSigmaAIR(Elev)

                # Compute Sigma UERE by combining all Sigma contributions
                # Ref: MOPS-DO-229D Section J.1
                #-----------------------------------------------------------------------
                # SCHEMATIC

                # Corrected Measurements from previous information
                #-----------------------------------------------------------------------
                # CorrPsr = SmoothedL1 + SatClk - UISD - STD

                # Compute the Geometrical Range
                # -----------------------------------------------------------------------
                # GeomRange = computeGeomRange(SATxyz, RCVRXYZ)

                # Compute the first Residual removing the geometrical range
                # -----------------------------------------------------------------------
                # PsrResidual = CorrPsr - GeomRange
            # End of if int(SatInfo[SatLabel][SatIdx['UDREI']]) < 12:

            # Prepare output for the satellite
            CorrInfo[SatLabel] = SatCorrInfo


        # End of if(SatPrepro["Status"] == 1):
    # End of for SatLabel, SatPrepro in PreproObsInfo.items():

    # Estimate the Receiver Clock first guess as a weighted average of the Residuals # (with the weights W=1/UERE2)
    # RcvrClk = estimateRcvrClk(PsrResidual, SigmaUERE)

    # Loop over satellites at PsrResidual Output
    # for Prn in Sat[Prn].Monitored:
        # Correct Residual from the first guess of the Receiver Clock
        # PsrResidual = PsrResidual - RcvrClk

    # Estimate the ENT-GPS Offset as the instantaneous average of the orbit corrections
    # #(LTC) projected into the user line-of-sight (Ulos) minus the clock corrections (FC+LTC)
    # ENTtoGPS = estimateENTtoGPS(FLT_Corrections)

    # SCHEMATIC

    return CorrInfo

# END OF runCorrectMeas()