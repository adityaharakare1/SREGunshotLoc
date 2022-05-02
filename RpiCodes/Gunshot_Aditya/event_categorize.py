import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def Friedlander(t,tRise,tDecay,tTrough):
    tPeak=0.
    Amp=1.
    PF = np.zeros_like(t)
    tCntr = t - tPeak
    itRise = np.where((-tRise <= tCntr) & (tCntr < 0))[0]
    if len(itRise):
        PF[itRise] = Amp*(1.+tCntr[itRise]/tRise)
    itDecay = np.where(tCntr >= 0)[0]
    PF[itDecay] = Amp*(1 - tCntr[itDecay]/tDecay) \
        *np.exp(-tCntr[itDecay]/(tTrough-tDecay))
    if len(t) == 1:
        PF = PF[0]
    return PF


def NWave(t,tRise,tDecay):
    tPeak=0.
    Amp=1.
    PS = np.zeros_like(t)
    tCntr = t - tPeak
    itRise1 = np.where((-tRise <= tCntr) & (tCntr < 0))[0]
    PS[itRise1] = Amp*(1.+tCntr[itRise1]/tRise)
    itDecay = np.where((0 <= tCntr) & (tCntr < tDecay))[0]
    PS[itDecay] = Amp*(1 - 2*tCntr[itDecay]/tDecay)
    itRise2 = np.where((tDecay <= tCntr) & (tCntr < tDecay+tRise))[0]
    PS[itRise2] = Amp*(tCntr[itRise2]-tDecay-tRise)/tRise
    return PS


def categorize(data, Fs, DEBUG):
    dt = 1/Fs
    _SW_Min_Neg_Peak = 0.5
    errors = np.zeros((4,2))
    for i in range(4):
        imax = np.argmax(data[i,:])  # Index of +ve peak of signal
        # print(imax)
        # Determine left and right end indices (relative to 'imax') to
        # window data (already windowed) based on '_TBefr' & '_TAftr'
        # ileft = max([imax - int(math.ceil(_TBefr/dt)),0])
        ileft = 0
        # iright = min([imax + int(math.ceil(_TAftr/dt)),len(snsr.data)])
        iright = len(data[i,:])
        # data = snsr.data[ileft:iright]  # Extract chosen window of data
        data[i,:] /= np.max(data[i,:])  # Normalize by max
        # Create time axis such that peak occurs at t=0
        tAxis = (np.arange(iright - ileft) - (imax - ileft)) * dt
        fitErrF = False  # As yet, no error encountered in Friedlander fit
        fitErrN = False  # As yet, no error encountered in N-wave fit
        # Generate initial guesses for various fit parameters
        jmax = np.argmax(data[i,:])  # Index of +ve peak of 'data'
        # Guess time to first 0-crossing working backward from +ve peak,
        # which is the rise time of the signal
        jzeroCrossMinus = np.argmax(data[i,:][jmax::-1] < 0)
        if jzeroCrossMinus == 0:  # No 0-crossing found in above
            t_Rise0 = -tAxis[0]  # Guess rise time is -ve of left time end
        else:
            t_Rise0 = -tAxis[jmax - jzeroCrossMinus]  # Estimate rise time
        # print(t_Rise0)
        try:  # Fitting Friedlander profile
            # Guess the time to first 0-crossing after +ve peak, which
            # is the decay time of the signal
            jzeroCrossPlus = np.argmax(data[i,:][jmax:] < 0)
            if jzeroCrossPlus == 0:  # No 0-crossings found in above
                MB_t_Decay0 = tAxis[-1]
            else:
                MB_t_Decay0 = tAxis[jmax + jzeroCrossPlus]
            # print(MB_t_Decay0)
            # Guess time from the +ve peak to the -ve peak of the signal
            jmin = np.argmin(data[i,:][jmax:])
            MB_t_Trough0 = tAxis[jmax + jmin]
            # print(MB_t_Trough0)
            if MB_t_Trough0 <= MB_t_Decay0:  # Can happen if no -ve swing
                MB_t_Trough0 = MB_t_Decay0 * 2  # Guess when -ve peak occurs
            # Curve fitting of sensor data using Friedlander function
            optF, covF = curve_fit(Friedlander, tAxis, data[i,:], p0=[t_Rise0, MB_t_Decay0, MB_t_Trough0])
            # print(optF)
        except:
            fitErrF = True  # Failure to fit
        else:
            # Successful fitting: calculate mean square error of fit
            FrFitData = Friedlander(tAxis, optF[0], optF[1], optF[2])
            if DEBUG:
                # Visualize Curve Fitting
                plt.plot(data[i,:])
                plt.plot(FrFitData, '-')
                plt.show()
            errors[i, 0] = np.average(np.square(data[i,:] - FrFitData))
            # print(errors[i,0])
        if np.min(data[i,:][jmax:]) < -_SW_Min_Neg_Peak:  # N-wave fit attempt?
            # (Normalized) -ve peak should be of sufficient magnitude to
            # warrant an attempt to fit an N-wave profile to the signal
            try:  # Fitting N-wave profile
                # Guess the time from the +ve peak to the -ve peak
                jmin = np.argmin(data[i,:])
                SW_t_Decay0 = tAxis[jmin]
                # Curve fitting of sensor data using N-wave function
                optN, covN = curve_fit(NWave, tAxis, data[i,:], p0=[t_Rise0, SW_t_Decay0])
            except:
                fitErrN = True  # Failure to fit
            else:
                # Successful fitting: calculate mean square error of fit
                NWFitData = NWave(tAxis, optN[0], optN[1])
                errors[i, 1] = np.average(np.square(data[i,:] \
                                                             - NWFitData))
        else:
            fitErrN = True  # Didn't even attempt an N-wave fit; => error
            errors[i, 1] = float('inf')


    minF = np.min(errors[:, 0])  # min error in Friedlander fits
    minN = np.min(errors[:, 1])  # min error in N-wave fits

    if minF != float('inf') or minN != float('inf'):  # At least one is finite
        if minF < minN:  # Friedlander fit was best across all sensors
            category = 'MuzzleBlast'  # Assign category
        else:  # N-wave fit was best across all sensors
            category = 'ShockWave' # Assign category
    else:
        # If both minima are inf (i.e. foregoing fits are all erroneous),
        # then category should be default
        category = -1
    print('Detected Event: ', category)
    return category
