import matplotlib.pyplot as plt
import numpy as np
import math

my_eps = 0.003


def calcDelayPair(data1, data2, Fs):
    """
    Calculates the DTOA between data 1 and data 2
    Args:
        data1: Data from microphone 1
        data2: Data from microphone 2
        Fs: Sampling Frequency

    Returns: Delay in secs

    """
    dt = 1/Fs
    corrArray = np.correlate(data1, data2, 'full')

    JDelay = np.argmax(corrArray)

    IDelay = (len(data1) - 1) - JDelay

    if JDelay == 0:
        corrPkNghbr = corrArray[:3]
    elif JDelay == len(corrArray) - 1:
        corrPkNghbr = corrArray[-3:]
    else:
        corrPkNghbr = corrArray[(JDelay - 1):(JDelay + 2)]
    if abs(corrPkNghbr[0] - 2 * corrPkNghbr[1] + corrPkNghbr[2]) > my_eps:
        DelFrctn = 0.5 * (corrPkNghbr[0] - corrPkNghbr[2]) \
                   / (corrPkNghbr[0] - 2 * corrPkNghbr[1] + corrPkNghbr[2])
    else:
        DelFrctn = 0.

    return (IDelay - DelFrctn) * dt


def calcDelaysAll(data, Fs):
    """
    Generates the DTOA Matrix
    Args:
        data: Data from 4 microphones
        Fs: Sampling Frequency

    Returns: Delays Matrux

    """
    Delays = np.zeros((4, 4), dtype=np.float)
    for iRef in range(4):
        for iWrk in range(iRef+1,4):
            Delays[iRef,iWrk] = calcDelayPair(data[iRef],data[iWrk],Fs)
            Delays[iWrk,iRef] = - Delays[iRef,iWrk]
    return Delays


def calSrcLoc(c0, c1, c2, Dist, Delays):
    sndSpd = 331.29
    d10 = Dist[c0, c1]
    d20 = Dist[c0, c2]
    # print (Delays[c0, c1]*80000)
    # print (Delays[c0, c2] * 80000)
    # normalized time taken by event to reach sensor 1 minus that to 0
    eta10 = Delays[c0, c1] * sndSpd / d10
    # normalized time taken by event to reach sensor 2 minus that to 0
    eta20 = Delays[c0, c2] * sndSpd / d20
    # print ('Normalized TDOA')
    # print (eta10, eta20)
    # if abs(eta10) >= 1.+my_eps or abs(eta20) >= 1.+my_eps:
    #     print ("Exception: Normalized TDOA > 1")
    #     return None

    thetaSensor10 = math.pi
    # determine angle between Sensors 0-1 arm and Sensors 0-2 arm
    Gamma = math.pi/4

    # print(my_eps)
    # find bearing angle (denoted 'srcTheta') w.r.t. Sensor 0-1 arm
    print(abs(eta10))
    if abs(eta10) > 1. - my_eps:  # source is collinear w/ Sensors 0-1
        print ('Src is colliner with Microphone 0-1')
        if eta10 > 0.:  # sound takes longer to reach Sensor 1 than Sensor 0
            # source and Sensor 1 are on opposite sides of Sensor 0
            srcThetas = [math.pi]
        else:  # sound takes longer to reach Sensor 0 than Sensor 1
            # soure and Sensor 1 are on same side of Sensor 0
            srcThetas = [0.]
    elif abs(eta20) > 1. - my_eps:  # source is collinear w/ Sensors 0-2
        print ('Src is colliner with Microphone 0-2')
        if eta20 > 0.:  # sound takes longer to reach Sensor 2 than Sensor 0
            # source and Sensor 2 are on opposite sides of Sensor 0
            srcThetas = [Gamma - math.pi]
        else:  # sound takes longer to reach Sensor 0 than Sensor 2
            # soure and Sensor 2 are on same side of Sensor 0
            srcThetas = [Gamma]
    else:  # source isn't collinear w/ either Sensors 0-1 or Sensors 0-2
        fctrA = (d20 / d10) * (1 - (eta20 ** 2)) / (1 - (eta10 ** 2))  # a factor
        alpha = math.atan2(math.sin(Gamma), math.cos(Gamma) - fctrA)  # an angle
        cosBeta = (fctrA * eta10 - eta20) / math.sqrt(1 + (fctrA ** 2) - (2 * fctrA * math.cos(Gamma)))  # cos of another angle
        if abs(cosBeta) > 1.:
            return None, None, None
        #            cosBeta = max(min(cosBeta,+1.),-1.)
        beta = math.acos(cosBeta)
        srcThetas = [alpha + beta, alpha - beta]  # two solutions

        # srcThetas = [3*math.pi/2 - (alpha + beta), 3*math.pi/2 -(alpha - beta)]

    print ("Possible Angles")
    for i in srcThetas:
        print (i*180/math.pi)

    # print([srcThetas[0]*180./math.pi, srcThetas[1]*180./math.pi])
    # internal function to get valid source location from bearing angle
    def getValidSoln(srcTheta):
        # use range expression from Sensors 0-1
        srcRange = d10 * (1 - eta10 ** 2) / (2 * (eta10 + math.cos(srcTheta)))
        # print(eta10 + math.cos(srcTheta))
        # print ("Range")
        # print (srcRange)
        # if srcRange < 0:
        #     return None, None
        # calculate absolute source angle
        srcThetaAbs = srcTheta + (thetaSensor10 - math.pi)
        # calculate absolute x-location of source
        srcLocX = srcRange * math.cos(srcThetaAbs)
        # calculate absolute y-location of source
        srcLocY = srcRange * math.sin(srcThetaAbs)
        # source's absolute z-location is preset as that of cluster centroid
        srcLocZ = 0
        srcLoc = np.array([srcLocX, srcLocY, srcLocZ])
        return srcThetaAbs, srcRange

    srcThetaAbs = []
    maxRange = -9999
    for srcTheta in srcThetas: #go thru all candidate solutions of srcTheta
        # valid source location solution from srcTheta ('None' if invalid)
        sol, range = getValidSoln(srcTheta)
        if sol is not None:
            if range>maxRange:
                srcThetaAbs.append(sol*180/math.pi)
                maxRange = range

    return srcThetaAbs[-1]

def calSrcLocAvg(data):
    r2 = math.sqrt(2)
    Fs = 80000
    # coords = {
    #     0: '0,0,0',
    #     1: '0,1,0',
    #     2: '1,0,0',
    #     3: '1,1,0'
    # }
    Dist = np.array([[0, 1, r2, 1], [1, 0, 1, r2], [r2, 1, 0, 1], [1, r2, 1, 0]])


    Delays = calcDelaysAll(data, Fs)
    # print(Delays)
    # print (Delays[0,1]*Fs) #t0-t1
    # print (Delays[0, 2] * Fs)
    # print (Delays[0, 3] * Fs)
    loc1 = calSrcLoc(0, 1, 2, Dist, Delays)
    print ('Angle of SRC = ', loc1)

