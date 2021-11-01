import sys, os

# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)

from collections import OrderedDict
from InputOutput import RcvrIdx
from COMMON import GnssConstants
from COMMON.Wlsq import wlsqComputation
import numpy as np

# The most relevant function here
def computeSpvtSolution(Conf, RcvrInfo, CorrInfo):
    """
    Objective: Compute the SBAS PVT Solution.

    Methodology:
        - Build the Geometry matrix G in XYZ and ENU
        - Build D as DOP-Matrix and PDOP
        - Build the weighting matrix W
        - Build the sparse matrix S
        - ...

    Assumptions:
        - It is assumed that the satellite position remains constant throughout the iterative process
        - It is assumed that the corrections (Iono, Tropo...) remain constant throughout the iterative process
        - The first receiver clock guess is set to zero
        - The first receiver position matches the exact receiver position

    :param Conf: Config info
    :param RcvrInfo: Receiver info
    :param CorrInfo: Corrected Measurements for current epoch
    :return: PosInfo, position info for current epoch
    """

    PosInfo = OrderedDict({})

    # Initialize some variables
    NumSatSol = 0
    GMatrixRows = []
    WMatrixDiagonal = []

    # Check if corrected information is available at current epoch
    if len(CorrInfo) > 0:
        PosInfo = {
            "Sod": 0.0,             # Second of day
            "Doy": 0,               # Day of year
            "Lon": 0.0,             # Receiver estimated longitude
            "Lat": 0.0,             # Receiver estimated latitude
            "Alt": 0.0,             # Receiver estimated altitude
            "Clk": 0.0,             # Receiver estimated clock
            "Sol": 0,               # 0: No solution 1: PA Sol 2: NPA Sol
            "NumSatVis": 0.0,       # Number of visible satellites
            "NumSatSol": 0.0,       # Number of visible satellites in solution
            "Hpe": 0.0,             # HPE
            "Vpe": 0.0,             # VPE
            "Epe": 0.0,             # EPE
            "Npe": 0.0,             # NPE
            "Hpl": 0.0,             # HPL
            "Vpl": 0.0,             # VPL
            "Hsi": 0.0,             # Horiontal Safety Index
            "Vsi": 0.0,             # Vertical Safety Index
            "Hdop": 0.0,            # HDOP
            "Vdop": 0.0,            # VDOP
            "Pdop": 0.0,            # PDOP
            "Tdop": 0.0,            # TDOP
        } # End of PosInfo

        # Get SoD and DoY
        PosInfo["Sod"] = CorrInfo[list(CorrInfo.keys())[0]]["Sod"]
        PosInfo["Doy"] = CorrInfo[list(CorrInfo.keys())[0]]["Doy"]
        # Get receiver coordinates first guess in the ENU reference frame
        PosInfo["Lon"] = float(RcvrInfo[RcvrIdx["LON"]])
        PosInfo["Lat"] = float(RcvrInfo[RcvrIdx["LAT"]])
        PosInfo["Alt"] = float(RcvrInfo[RcvrIdx["ALT"]])

        # Traverse all valid satellites within corrected measurements
        for SatLabel, SatCorrInfo in CorrInfo.items():
            # PA - Precision Approach Service
            if SatCorrInfo["Flag"] == 1:
                # Update number of available satellites for PA
                NumSatSol = NumSatSol + 1
                # Build the G - GEOMETRY matrix
                GMatrixRows.append(buildGMatrix(SatCorrInfo))
                # Build the W - WEIGHT matrix
                WMatrixDiagonal.append(buildWMatrix(SatCorrInfo))
        # Visible Satellites
        PosInfo["NumSatVis"] = len(CorrInfo)
        # Visible Satellites && Used in the Solution
        PosInfo["NumSatSol"] = NumSatSol

        # SBAS PVT Computation
        # ----------------------------------------------------------------------
        # Compute inputs required for the computation
        # Build the G - GEOMETRY matrix
        GMatrix = np.reshape(GMatrixRows, (NumSatSol, 4))
        # Build the W - WEIGHT matrix
        WMatrix = np.diag(WMatrixDiagonal)

        # Check the number of available satellites for computing the solution
        if NumSatSol >= GnssConstants.MIN_NUM_SATS_PVT:
            # Compute DOPS
            computeDops(GMatrix, PosInfo)
            if PosInfo["Pdop"] < float(Conf["PDOP_MAX"]):
                # Compute the S matrix
                SMatrix = buildSMatrix(GMatrix, WMatrix)
                # Call WLSQ function
                wlsqComputation(Conf, CorrInfo, PosInfo, SMatrix)
                # Compute protection levels
                computeProtectionLevels(GMatrix, WMatrix, PosInfo)
                # Compute safety indexes
                computeSafetyIndexes(PosInfo)

                # Set spvt computed as PA solution ---- Better do it within the wlsq solver
                # PosInfo['Sol'] = 1

            else:
                # No spvt solution
                PosInfo['Sol'] = 0
        else:
            # No spvt solution
            PosInfo['Sol'] = 0

    return PosInfo


def computeSafetyIndexes(PosInfo):
    PosInfo["Hsi"] = PosInfo["Hpe"] / PosInfo["Hpl"]
    PosInfo["Vsi"] = PosInfo["Vpe"] / PosInfo["Vpl"]


def buildWMatrix(SatCorrInfo):
    """
    Compute diagonal element of the weighting matrix W for the spvt computation
    :param SatCorrInfo:
    :return:
    """

    # Compute the diagonal element of the weighting matrix W
    Weight = 1/(SatCorrInfo["SigmaUere"])**2

    return Weight

def buildGMatrix(SatCorrInfo):
    """
    Compute row of the geometry matrix G in the ENU reference frame for the Spvt computation
    :param SatCorrInfo:
    :return:
    """

    # Convert degrees to radians
    DegtoRad = np.pi / 180.0

    # Compute row of geometry matrix G
    XComp = - (np.cos(SatCorrInfo["Elevation"] * DegtoRad) * np.sin(SatCorrInfo["Azimuth"] * DegtoRad))
    YComp = - (np.cos(SatCorrInfo["Elevation"] * DegtoRad) * np.cos(SatCorrInfo["Azimuth"] * DegtoRad))
    ZComp = - (np.sin(SatCorrInfo["Elevation"] * DegtoRad))
    BComp = 1

    return XComp, YComp, ZComp, BComp

def buildSMatrix(GMatrix, WMatrix):
    """
    Compute the S matrix in the ENU reference frame for the svpt computation.
    """
    return np.linalg.multi_dot(
        [np.linalg.inv(
            np.linalg.multi_dot([GMatrix.T, WMatrix, GMatrix])), GMatrix.T, WMatrix])

def computeDops(GMatrix, PosInfo):
    """
    Computes the DOP matrix Q in the ENU reference frame for the spvt computation.
    Also computes HDOP, VDOP, PDOP and TDOP.

    :param GMatrix: 
    :param PosInfo: 
    :return: 
    """
    # Compute the DOP matrix
    QMatrix = np.linalg.inv(np.dot(GMatrix.T, GMatrix))
    QDiag = np.diag(QMatrix)

    # Compute the DOPS
    PosInfo["Hdop"] = np.sqrt(QDiag[0] + QDiag[1])
    PosInfo["Vdop"] = np.sqrt(QDiag[2])
    PosInfo["Pdop"] = np.sqrt(QDiag[0] + QDiag[1] + QDiag[2])
    PosInfo["Tdop"] = np.sqrt(QDiag[3])

def computeProtectionLevels(GMatrix, WMatrix, PosInfo):
    """
    Computes the EGNOS protection levels
    Ref Sys: ENU
    :param GMatrix:
    :param WMatrix:
    :param PosInfo:
    :return:
    """

    # D Matrix
    DMatrix = np.linalg.inv(np.linalg.multi_dot([GMatrix.T, WMatrix, GMatrix]))
    DMatrix_diag1 = np.diag(DMatrix)
    DMatrix_diag2 = np.diag(DMatrix, k = 1)

    # Compute the protection levels (PA service)
    PosInfo["Hpl"] = np.sqrt(((DMatrix_diag1[0] + DMatrix_diag1[1])/2) + np.sqrt(((DMatrix_diag1[0] - DMatrix_diag1[1])/2)**2 + DMatrix_diag2[0]**2))*GnssConstants.MOPS_KH_PA
    PosInfo["Vpl"] = np.sqrt(DMatrix_diag1[2])*GnssConstants.MOPS_KV_PA
