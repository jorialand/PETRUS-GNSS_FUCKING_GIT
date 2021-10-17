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

def buildResidualsVector(CorrInfo, PosInfo, Mode):

    # Purpose: Compute the residuals vector per iteration for the wlsq computation

    # Initialize the residuals vector 
    ResVector = []

    # Get receiver position in the WGS84 reference frame
    RcvrPos = np.array(llh2xyz(PosInfo["Lon"], PosInfo["Lat"], PosInfo["Alt"]))

    # Loop over all satellites in CorrInfo
    for SatCorrInfo in CorrInfo.values():
        if SatCorrInfo["Flag"] == 1:
            # Obtain satellite position in the WGS84 reference frame
            SatPos = np.array([SatCorrInfo["SatX"], SatCorrInfo["SatY"], SatCorrInfo["SatZ"]])
            # Compute the geometrical range
            GeomRange = np.linalg.norm(np.subtract(SatPos, RcvrPos))
            # Compute the residuals
            ResVector.append(SatCorrInfo["CorrPsr"] - PosInfo["Clk"] - GeomRange)
        if SatCorrInfo["Flag"] == 2 and Mode == "NPA":
            # Obtain satellite position in the WGS84 reference frame
            SatPos = np.array([SatCorrInfo["SatX"], SatCorrInfo["SatY"], SatCorrInfo["SatZ"]])
            # Compute the geometrical range
            GeomRange = np.linalg.norm(np.subtract(SatPos, RcvrPos))
            # Compute the residuals
            ResVector.append(SatCorrInfo["CorrPsr"] - PosInfo["Clk"] - GeomRange)

    # Get full residuals vector for the svpt computation
    ResVector = np.reshape(ResVector,(PosInfo["NumSatSol"],1))

    return ResVector

# Wlsq iterative computation
#-----------------------------------------------------------------------

def wlsqComputation(Conf, CorrInfo, PosInfo, SMatrix, Mode):

    # Purpose: Perform the wlsq iterative process as part of the svpt computation

    # More in detail, this function handles the following tasks:

    #       *  Build the residuals vector for each iteration
    #       *  Compute the final estimated position of the receiver
    #       *  Compute the position errors

    # Parameters
    # ==========
    # Conf: dict
    #       Configuration dictionary
    # CorrInfo: dict
    #           Corrected information for current epoch per satellite
    #           CorrInfo["G01"]["C1"]
    # PosInfo: dict
    #          Spvt outputs for Rcvr for current epoch: Latitude, longitude...
    # SMatrix: np.array
    #          Matrix used for the wlsq computation
    #          
    
    # Returns
    # =======
    # Nothing

    # Initialize some variables
    NumIter = 0
    NormPosDelta = 1000.0
    DegtoRad = np.pi/180.0

    # Wlsq computation
    while NumIter <= Conf["MAX_LSQ_ITER"] and NormPosDelta > GnssConstants.LSQ_DELTA_EPS:
        # Increase the number of iterations
        NumIter = NumIter + 1
        # Compute the residuals vector
        ResVector = buildResidualsVector(CorrInfo, PosInfo, Mode)
        # Compute the position residuals in the ENU reference frame 
        PosDelta = np.dot(SMatrix, ResVector).flatten()
        NormPosDelta = np.linalg.norm(PosDelta)
        # Update receiver position in geodesic coordinates 
        PosInfo["Lat"] = PosInfo["Lat"] + PosDelta[1]/(GnssConstants.EARTH_RADIUS*DegtoRad)
        PosInfo["Lon"] = PosInfo["Lon"] + PosDelta[0]/(GnssConstants.EARTH_RADIUS*np.cos(PosInfo["Lat"]*DegtoRad)*DegtoRad)
        PosInfo["Alt"] = PosInfo["Alt"] + PosDelta[2]
        PosInfo["Clk"] = PosInfo["Clk"] + PosDelta[3]
        # Estimate ENU position errors
        PosInfo["Epe"] = PosInfo["Epe"] + PosDelta[0]
        PosInfo["Npe"] = PosInfo["Npe"] + PosDelta[1]
        PosInfo["Hpe"] = np.sqrt(PosInfo["Epe"]**2 + PosInfo["Npe"]**2)
        PosInfo["Vpe"] = PosInfo["Vpe"] + PosDelta[2]

    # End of while(NumIter <= Conf["MAX_LSQ_ITER"] and NormPosDelta > GnssConstants.LSQ_DELTA_EPS):

    # If the wlsq iterative filter converges
    if NumIter <= Conf["MAX_LSQ_ITER"]:
        # PA solution achieved
        if Mode == "PA":
            PosInfo["Sol"] = 1
        # NPA solution achieved
        elif Mode == "NPA":
            PosInfo["Sol"] = 2

    # End of if(NumIter <= Conf["MAX_LSQ_ITER"]):

# End of wlsq computation
# -----------------------------------------------------------------------