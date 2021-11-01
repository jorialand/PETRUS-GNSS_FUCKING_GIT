import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from COMMON import GnssConstants
from COMMON.Coordinates import llh2xyz
import numpy as np

# Wlsq internal functions
#-----------------------------------------------------------------------

def BuildPseudorangeResiduals(CorrInfo, PosInfo):
    """
    Update Pseudorange Residuals
    :param CorrInfo:
    :param PosInfo:
    :param Mode:
    :return:
    """

    PsrResiduals = []
    # Rcvr position (Ref Sys: WGS84)
    RcvrPos = np.array(llh2xyz(PosInfo["Lon"], PosInfo["Lat"], PosInfo["Alt"]))

    for SatCorrInfo in CorrInfo.values():
        if SatCorrInfo["Flag"] == 1:
            # SV position (Ref Sys: WGS84)
            SatPos = np.array([SatCorrInfo["SatX"], SatCorrInfo["SatY"], SatCorrInfo["SatZ"]])
            GeomRange = np.linalg.norm(np.subtract(SatPos, RcvrPos))

            # Compute the residuals
            PsrResiduals.append(SatCorrInfo["CorrPsr"] - PosInfo["Clk"] - GeomRange)
        else:
            # Non converged measurement, or service level different than PA
            # sys.stderr.write("WARNING: Expected  PA measurement, got something else instead.\n")
            pass

    # Get full residuals vector for the svpt computation
    PsrResiduals = np.reshape(PsrResiduals,(PosInfo["NumSatSol"],1))

    return PsrResiduals

# Wlsq iterative computation
#-----------------------------------------------------------------------

def wlsqComputation(Conf, CorrInfo, PosInfo, SMatrix):
    """
    Weighted Least Squares solver
    """

    # Initialize some variables
    iter = 0
    RcvrPosDelta = 999

    # Wlsq computation
    while iter < Conf["MAX_LSQ_ITER"] and RcvrPosDelta > GnssConstants.LSQ_DELTA_EPS:
        # Update Pseudorange Residuals
        PsrResiduals = BuildPseudorangeResiduals(CorrInfo, PosInfo)

        # Compute the position residuals in the ENU reference frame 
        PosDelta = np.dot(SMatrix, PsrResiduals).flatten()

        # Rcvr position residuals
        RcvrPosDelta = np.linalg.norm(PosDelta)
        # Update Rcvr pos in geodetic coords
        PosInfo["Lat"] = PosInfo["Lat"] + np.rad2deg(PosDelta[1]/(GnssConstants.EARTH_RADIUS))
        PosInfo["Lon"] = PosInfo["Lon"] + np.rad2deg(PosDelta[0]/(GnssConstants.EARTH_RADIUS*np.cos(np.deg2rad(PosInfo["Lat"]))))
        PosInfo["Alt"] = PosInfo["Alt"] + PosDelta[2]
        PosInfo["Clk"] = PosInfo["Clk"] + PosDelta[3]
        # Estimate ENU position errors
        PosInfo["Epe"] = PosInfo["Epe"] + PosDelta[0]
        PosInfo["Npe"] = PosInfo["Npe"] + PosDelta[1]
        PosInfo["Hpe"] = np.linalg.norm([PosInfo['Epe'], PosInfo['Npe']])
        PosInfo["Vpe"] = PosInfo["Vpe"] + PosDelta[2]

        # Increase the number of iterations
        iter += 1

    # End of while(NumIter <= Conf["MAX_LSQ_ITER"] and NormPosDelta > GnssConstants.LSQ_DELTA_EPS):

    # If the wlsq iterative filter converges
    if iter < Conf["MAX_LSQ_ITER"]:
        # PA solution achieved
        PosInfo["Sol"] = 1