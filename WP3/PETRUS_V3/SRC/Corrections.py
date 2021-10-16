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

IgpIdx2Vertex = {
    1: "NE",
    2: "NW",
    3: "SW",
    4: "SE",
}

IgpVertex2Idx = {
    "NE": 1,
    "NW": 2,
    "SW": 3,
    "SE": 4,
}

def correctSatPosAndClk(SatInfo, CorrectInfo):
    # Reference: MOPS-DO-229D Section A.4.4.7

    # Apply LTC-X
    CorrectInfo["SatX"] = float(SatInfo[SatIdx["SAT-X"]]) + float(SatInfo[SatIdx["LTC-X"]])
    
    # Apply LTC-Y
    CorrectInfo["SatY"] = float(SatInfo[SatIdx["SAT-Y"]]) + float(SatInfo[SatIdx["LTC-Y"]])
    
    # Apply LTC-Z
    CorrectInfo["SatZ"] = float(SatInfo[SatIdx["SAT-Z"]]) + float(SatInfo[SatIdx["LTC-Z"]])
    
    # Compute Dtr
    Dtr = (-2 * \
            np.dot(\
                [float(SatInfo[SatIdx["SAT-X"]]),
                float(SatInfo[SatIdx["SAT-Y"]]),
                float(SatInfo[SatIdx["SAT-Z"]])],
                [float(SatInfo[SatIdx["VEL-X"]]),
                float(SatInfo[SatIdx["VEL-Y"]]),
                float(SatInfo[SatIdx["VEL-Z"]])],
            )\
        ) / Const.SPEED_OF_LIGHT

    # Apply clock corrections
    CorrectInfo["SatClk"] = float(SatInfo[SatIdx["SAT-CLK"]]) + \
        (-1)*float(SatInfo[SatIdx["TGD"]]) + \
        Dtr + \
        float(SatInfo[SatIdx["FC"]]) + \
        float(SatInfo[SatIdx["LTC-B"]]) \


def computeSigmaFlt(SatInfo, CorrectInfo):
    # Reference: MOPS-DO-229D Section A.4.5.1

    # If RSS==0
    if(SatInfo[SatIdx["RSS"]]=='0'):
        # Compute non-root-sum-squared
        CorrectInfo["SigmaFlt"] = (\
            (float(SatInfo[SatIdx["SIGMAUDRE"]]) * float(SatInfo[SatIdx["DELTAUDRE"]])) + \
                float(SatInfo[SatIdx["EPS-FC"]]) + \
                float(SatInfo[SatIdx["EPS-RRC"]]) + \
                float(SatInfo[SatIdx["EPS-LTC"]]) + \
                float(SatInfo[SatIdx["EPS-ER"]])
                    )

    else:
        # Compute root-sum-squared
        CorrectInfo["SigmaFlt"] = np.sqrt(\
            (float(SatInfo[SatIdx["SIGMAUDRE"]]) * float(SatInfo[SatIdx["DELTAUDRE"]]))**2 +\
                (float(SatInfo[SatIdx["EPS-FC"]]))**2 + \
                (float(SatInfo[SatIdx["EPS-RRC"]]))**2 + \
                (float(SatInfo[SatIdx["EPS-LTC"]]))**2 + \
                (float(SatInfo[SatIdx["EPS-ER"]]))**2
        )


def rewrapLon(Longitude):
    if abs(Longitude) > 180.0:
        return (Longitude - \
                (Longitude/abs(Longitude)) * 360.0
                )

    else:
        return (Longitude)


def rewrapLat(Latitude):
    if abs(Latitude) > 90.0:
        return (Latitude - \
                (Latitude/abs(Latitude)) * 180.0
                )

    else:
        return (Latitude)


def computeUisdAndUire(PreproInfo, LosInfo, CorrectInfo):
    # Reference: MOPS-DO-229D Section A.4.4.10.3

    # Rectangular interpolation
    if(LosInfo[LosIdx["INTERP"]] == '0'):
        # Case of IPP between S85 and N85
        if (float(LosInfo[LosIdx["IPPLAT"]]) < 85.0) and\
            (float(LosInfo[LosIdx["IPPLAT"]]) > -85.0):
            # Compute xpp
            xpp = rewrapLon(float(LosInfo[LosIdx["IPPLON"]]) - float(LosInfo[LosIdx["IGP_SW_LON"]]))/\
                rewrapLon(float(LosInfo[LosIdx["IGP_SE_LON"]]) - float(LosInfo[LosIdx["IGP_SW_LON"]]))

            # Compute ypp
            ypp = rewrapLat(float(LosInfo[LosIdx["IPPLAT"]]) - float(LosInfo[LosIdx["IGP_SW_LAT"]]))/\
                rewrapLat(float(LosInfo[LosIdx["IGP_NW_LAT"]]) - float(LosInfo[LosIdx["IGP_SW_LAT"]]))

        else:
            # Compute ypp
            ypp = (abs(float(LosInfo[LosIdx["IPPLAT"]])) - 85.0) / 10.0
            
            # Compute Delta Longitude
            # Compute the intermediary coefficients
            DeltaLon = rewrapLon(\
                float(LosInfo[LosIdx["IPPLON"]]) - \
                float(LosInfo[LosIdx["IGP_SW_LON"]])
            )
            
            # Compute xpp
            xpp = ((DeltaLon/90.0) * (1.0 - (2.0 * ypp))) + ypp 

        # End of if (float(LosInfo[LosIdx["IPPLAT"]]) < 85.0) and\...

        # Compute the interpolation weights
        W1 = xpp * ypp
        W2 = (1 - xpp) * ypp
        W3 = (1 - xpp) * (1 - ypp)
        W4 = xpp * (1 - ypp)

        # Compute UISD
        CorrectInfo["Uisd"] = PreproInfo["Mpp"] * (\
            (W1 * float(LosInfo[LosIdx["GIVD_NE"]])) +\
            (W2 * float(LosInfo[LosIdx["GIVD_NW"]])) +\
            (W3 * float(LosInfo[LosIdx["GIVD_SW"]])) +\
            (W4 * float(LosInfo[LosIdx["GIVD_SE"]]))
        )

        # Compute UIRE
        CorrectInfo["SigmaUire"] = np.sqrt(\
            PreproInfo["Mpp"]**2 * (\
                (W1 * float(LosInfo[LosIdx["GIVE_NE"]])**2) +\
                (W2 * float(LosInfo[LosIdx["GIVE_NW"]])**2) +\
                (W3 * float(LosInfo[LosIdx["GIVE_SW"]])**2) +\
                (W4 * float(LosInfo[LosIdx["GIVE_SE"]])**2)
            )
        )

    # Triangular interpolation
    else:
        # Get index of the Vertex opposite the hypotenuse
        Idx2 = (int(LosInfo[LosIdx["INTERP"]]) + 2) % 4
        if Idx2==0: Idx2=4
        # Get Vertex opposite the hypotenuse
        Vertex2 = IgpIdx2Vertex[Idx2]

        # Get Vertex 1
        Vertex1 = "%s" % Vertex2
        Vertex1 = (Vertex1.replace('S', 'N') if 'S' in Vertex1 else Vertex1.replace('N', 'S'))

        # Get Vertex 3
        Vertex3 = "%s" % Vertex2
        Vertex3 = (Vertex3.replace('E', 'W') if 'E' in Vertex3 else Vertex3.replace('W', 'E'))

        # Compute xpp
        xpp = rewrapLon(float(LosInfo[LosIdx["IPPLON"]]) - \
            float(LosInfo[LosIdx["IGP_" + Vertex2 + "_LON"]]))/\
            rewrapLon(float(LosInfo[LosIdx["IGP_" + Vertex3 + "_LON"]]) - \
                float(LosInfo[LosIdx["IGP_" + Vertex2 + "_LON"]]))

        # Compute ypp
        ypp = rewrapLat(float(LosInfo[LosIdx["IPPLAT"]]) - \
            float(LosInfo[LosIdx["IGP_" + Vertex2 + "_LAT"]]))/\
                rewrapLat(float(LosInfo[LosIdx["IGP_" + Vertex1 + "_LAT"]]) - \
                    float(LosInfo[LosIdx["IGP_" + Vertex2 + "_LAT"]]))

        # Compute the interpolation weights
        W1 = ypp
        W2 = 1 - xpp - ypp
        W3 = xpp

        # Compute UISD
        CorrectInfo["Uisd"] = PreproInfo["Mpp"] * (\
            (W1 * float(LosInfo[LosIdx["GIVD_" + Vertex1]])) +\
            (W2 * float(LosInfo[LosIdx["GIVD_" + Vertex2]])) +\
            (W3 * float(LosInfo[LosIdx["GIVD_" + Vertex3]]))
        )

        # Compute UIRE
        CorrectInfo["SigmaUire"] = np.sqrt(\
            PreproInfo["Mpp"]**2 * (\
                (W1 * float(LosInfo[LosIdx["GIVE_" + Vertex1]])**2) +\
                (W2 * float(LosInfo[LosIdx["GIVE_" + Vertex2]])**2) +\
                (W3 * float(LosInfo[LosIdx["GIVE_" + Vertex3]])**2)
            )
        )

    # End of if(int(LosInfo[LosIdx["INTERP"]]) == '0'):


def computeTropoMpp(Elev):
    if Elev>=4:
        TropoMpp = (1.001/(np.sqrt(0.002001+(np.sin(np.radians(Elev)))**2)))
    elif Elev>=2: 
        TropoMpp = (1.001/(np.sqrt(0.002001+(np.sin(np.radians(Elev)))**2)))*\
        (1+0.015*(max(0,4-Elev))**2)
    else:
        TropoMpp = np.nan

    return TropoMpp


def computeSigmaTropo(TropoMpp):
    SigmaTropo = 0.12*TropoMpp

    return SigmaTropo

def computeSigmaAirborne(Conf, Elev, CorrectInfo):
    if Conf["EQUIPMENT_CLASS"] ==1:
        SigmaAirborne = 5

    else:
        SigmaMpSquare = (0.13+0.53*np.exp(-Elev/10.0))**2
        if Conf["AIR_ACC_DESIG"] == 'A':
            if Elev > Conf["ELEV_NOISE_TH"]:
                SigmaNoiseDivSquare = 0.15 ** 2

            else:
                SigmaNoiseDivSquare = 0.36 ** 2

        elif Conf["AIR_ACC_DESIG"] == 'B':
            if Elev > Conf["ELEV_NOISE_TH"]:
                SigmaNoiseDivSquare = 0.11 ** 2

            else:
                SigmaNoiseDivSquare = 0.15 ** 2

        SigmaAirborne = np.sqrt(SigmaMpSquare + SigmaNoiseDivSquare)

    CorrectInfo["SigmaNoiseDiv"] = np.sqrt(SigmaNoiseDivSquare)
    CorrectInfo["SigmaMultipath"] = np.sqrt(SigmaMpSquare)
    CorrectInfo["SigmaAirborne"] = SigmaAirborne


def computeEntGps(SatInfo, Rcvr):
    LtcXyz =  np.array(\
        [float(SatInfo[SatIdx["LTC-X"]]),
        float(SatInfo[SatIdx["LTC-Y"]]),
        float(SatInfo[SatIdx["LTC-Z"]])])

    Ulos =  np.array(\
        [float(SatInfo[SatIdx["SAT-X"]]) - Rcvr[RcvrIdx["XYZ"]][0],
        float(SatInfo[SatIdx["SAT-Y"]]) - Rcvr[RcvrIdx["XYZ"]][1],
        float(SatInfo[SatIdx["SAT-Z"]]) - Rcvr[RcvrIdx["XYZ"]][2]])
    Unorm = np.linalg.norm(Ulos)
    Ulos = Ulos/Unorm

    EntGps = (np.dot(LtcXyz, Ulos) - \
        (float(SatInfo[SatIdx["FC"]]) + \
            float(SatInfo[SatIdx["LTC-B"]]))\
    )

    return EntGps

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
                "SatX": 0.0,            # X-Component of the Satellite Position 
                                        # corrected with SBAS LTC
                "SatY": 0.0,            # Y-Component of the Satellite Position 
                                        # corrected with SBAS LTC
                "SatZ": 0.0,            # Z-Component of the Satellite Position 
                                        # corrected with SBAS LTC
                "SatClk": 0.0,          # Satellite Clock corrected with SBAS FLT
                "Uisd": 0.0,            # User Ionospheric Slant Delay
                "Std": 0.0,             # Slant Tropospheric Delay
                "CorrPsr": 0.0,         # Pseudo Range corrected from delays
                "GeomRange": 0.0,       # Geometrical Range (distance between Satellite 
                                        # Position and Receiver Reference Position)
                "PsrResidual": 0.0,     # Pseudo Range Residual
                "RcvrClk": 0.0,         # Receiver Clock estimation
                "SigmaFlt": 0,          # Sigma of the residual error associated to the 
                                        # fast and long-term correction (FLT)
                "SigmaUire": 0,         # User Ionospheric Range Error Sigma
                "SigmaTropo": 0,        # Sigma of the Tropospheric error 
                "SigmaAirborne": 0.0,   # Sigma Airborne Error
                "SigmaNoiseDiv": 0.0,   # Sigma of the receiver noise + divergence
                "SigmaMultipath": 0.0,  # Sigma of the receiver multipath
                "SigmaUere": 0.0,       # Sigma User Equivalent Range Error (Sigma of 
                                        # the total residual error associated to the 
                                        # satellite)
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
            if (SatLabel in SatInfo) and (SatLabel in LosInfo):
                # Get IPP Longitude
                SatCorrInfo["IppLon"] = float(LosInfo[SatLabel][LosIdx["IPPLON"]])
                # Get IPP Latitude
                SatCorrInfo["IppLat"] = float(LosInfo[SatLabel][LosIdx["IPPLAT"]])

                # If satellite is Not Monitored or Don't Use, continue to next satellite
                if(int(SatInfo[SatLabel][SatIdx["UDREI"]]) >= 14):
                    # Set LoS flag to 0
                    SatCorrInfo["Flag"] = 0

                    # Prepare output for the satellite
                    CorrInfo[SatLabel] = SatCorrInfo

                    continue

                elif(int(SatInfo[SatLabel][SatIdx["UDREI"]]) >= 12):
                    # Set LoS flag to NPA
                    SatCorrInfo["Flag"] = 2

                # End of if(int(SatInfo[SatLabel][SatIdx["UDREI"]]) >= 14):

                # Apply the SBAS corrections to the satellite position and clock
                correctSatPosAndClk(SatInfo[SatLabel], SatCorrInfo)

                # Compute the Sigma FLT projected into the User direction as per MOPS
                computeSigmaFlt(SatInfo[SatLabel], SatCorrInfo)

                # # Compute UISD and UIRE on the IPP using MOPS interpolation (Appendix A) (TODO)
                computeUisdAndUire(SatPrepro, LosInfo[SatLabel], SatCorrInfo)
                # SatCorrInfo["Uisd"] = float(LosInfo[SatLabel][LosIdx["UISD"]])
                # SatCorrInfo["SigmaUire"] = float(LosInfo[SatLabel][LosIdx["SUIRE"]])

                # Compute the STD: Slant Tropo Delay and associated SigmaTROPO
                # Refer to MOPS guidelines in Appendix A section A.4.2.4 for Tropospheric 
                # Model
                #-----------------------------------------------------------------------
                # Compute Tropospheric Mapping Function
                TropoMpp = computeTropoMpp(SatCorrInfo["Elevation"])

                # # [OPTIONAL] Compute the Slant Tropospheric Delay (TODO)
                # SatCorrInfo["Std"] = computeSlantTropoDelay(RCVR[iRec].llh, Doy)
                SatCorrInfo["Std"] = float(LosInfo[SatLabel][LosIdx["STD"]])
                
                # Compute the Slant Tropospheric Delay Error Sigma
                SatCorrInfo["SigmaTropo"] = computeSigmaTropo(TropoMpp)

                # Compute User Airborne Sigma
                # Ref: MOPS-DO-229D Section J.2.4
                #-----------------------------------------------------------------------
                # Compute SigmaAIR: Sigma Airborne in line with SBAS Standard
                computeSigmaAirborne(Conf, SatCorrInfo["Elevation"], SatCorrInfo)

                # Compute UERE by combining all Sigma contributions
                # Ref: MOPS-DO-229D Section J.1
                #-----------------------------------------------------------------------
                # Compute UERE
                SatCorrInfo["SigmaUere"] = np.sqrt(\
                        SatCorrInfo["SigmaFlt"]**2 + \
                        SatCorrInfo["SigmaUire"]**2 + \
                        SatCorrInfo["SigmaTropo"]**2 + \
                        SatCorrInfo["SigmaAirborne"]**2 \
                )

                # Correct the Smoothed Pseudo Range from Sat Clock, Tropo and Iono delays
                #-----------------------------------------------------------------------
                SatCorrInfo["CorrPsr"] = \
                    SatPrepro["SmoothC1"] + SatCorrInfo["SatClk"] - \
                        SatCorrInfo["Uisd"] - SatCorrInfo["Std"]

                # ###########################################
                # Theta = Const.OMEGA_EARTH*abs(0.07)
                # RotationMatrix = \
                #     [
                #         [np.cos(Theta), np.sin(Theta), 0],\
                #         [-np.sin(Theta), np.cos(Theta), 0],\
                #         [0,0,1]\
                #     ]

                # RotatedSatPos = \
                #     (np.dot(RotationMatrix, np.array([SatCorrInfo["SatX"], SatCorrInfo["SatY"], SatCorrInfo["SatZ"]])))

                # SatCorrInfo["SatX"] = RotatedSatPos[0]
                # SatCorrInfo["SatY"] = RotatedSatPos[1]
                # SatCorrInfo["SatZ"] = RotatedSatPos[2]
                # ###########################################

                # Compute the Geometrical Range
                SatCorrInfo["GeomRange"] = np.sqrt(\
                    (SatCorrInfo["SatX"] - Rcvr[RcvrIdx["XYZ"]][0])**2 +
                    (SatCorrInfo["SatY"] - Rcvr[RcvrIdx["XYZ"]][1])**2 +
                    (SatCorrInfo["SatZ"] - Rcvr[RcvrIdx["XYZ"]][2])**2
                )

                # Compute the Residual including Receiver Clock estimation
                SatCorrInfo["PsrResidual"] = \
                    SatCorrInfo["CorrPsr"] -  SatCorrInfo["GeomRange"]

                # Update the parameters to compute the Receiver Clock estimation
                ResSum = ResSum + ((SatCorrInfo["SigmaUere"]**-2) * SatCorrInfo["PsrResidual"])
                ResN = ResN + (SatCorrInfo["SigmaUere"]**-2)

                # Compute ENT-GPS estimation from current satellite
                EntGps = computeEntGps(SatInfo[SatLabel], Rcvr)

                # Update the parameters to compute the ENT-GPS Offset
                EntGpsSum = EntGpsSum + EntGps
                EntGpsN = EntGpsN + 1

            else:
                # Set LoS flag to 0
                SatCorrInfo["Flag"] = 0

            # End of if SatLabel in SatInfo

            # Prepare output for the satellite
            CorrInfo[SatLabel] = SatCorrInfo

        # End of if(SatPrepro["Status"] == 1):

    # End of for SatLabel, SatPrepro in PreproObsInfo.items():

    # Loop over corrected measurements
    for SatLabel, SatCorrInfo in CorrInfo.items():
        # Check if FLAG is set to 0
        if(SatCorrInfo["Flag"] > 0):
            # Compute the Receiver Clock estimation
            SatCorrInfo["RcvrClk"] = ResSum / ResN if ResN else np.nan

            # Correct Residuals from Receiver Clock estimation
            SatCorrInfo["PsrResidual"] = \
                SatCorrInfo["PsrResidual"] - SatCorrInfo["RcvrClk"]

            # Compute the ENT-GPS Offset
            SatCorrInfo["EntGps"] = EntGpsSum / EntGpsN if EntGpsN else np.nan

    return CorrInfo
