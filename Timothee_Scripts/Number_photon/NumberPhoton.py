# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 09:32:09 2020

@author: TimGU
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def RadToMHz(rad):
    return rad  / (2*np.pi)

def MHzToRad(mhz):
    return mhz * (2*np.pi)

def nbOsci(dist_vec, bg_dist):
    nbOsci = 0
    for i in range(3, len(dist_vec)):
        s0 = dist_vec[i-3] > bg_dist
        s1 = dist_vec[i-2] > bg_dist
        s2 = dist_vec[i-1] > bg_dist
        s3 = dist_vec[i] > bg_dist
        if s0==s1 and s2==s3 and s1!=s3:
            nbOsci +=1
    return nbOsci

def Periode(timeVec, signal, bg):
    peaks, _ = find_peaks(signal, width=5, height=bg)
    if len(peaks)>=2:
        return timeVec[peaks[1]] - timeVec[peaks[0]]
    
    nb_osci = nbOsci(signal, bg) 
    if nb_osci != 0:
        return timeVec[-1] / nbOsci(signal, bg) 
    
    return timeVec[-1] / 5
    

###############################################
######## FIT FUNCTIONS ########################
###############################################

def ramseyFit(t, bg, osci, Wr, phi0, Tr):
    return bg - osci*np.sin(Wr*t + phi0)*np.exp(-t/Tr)

def ramseyStarkDephasing(x, bg, wc, k, chi, ed):
    num = 8*k* chi**2 * ed
    den = ( k**2 + 4*(x-wc)**2 - chi**2 )**2   +  4* chi**2 * k**2
    
    return bg + num/den

def ramseyStarkQbShift(x, bg, wc, k, chi, ed):
    num = 4 * chi * ed * (  k**2 + 4*(x-wc)**2 - chi**2  )
    den = ( k**2 + 4*(x-wc)**2 - chi**2 )**2   +  4* chi**2 * k**2
    return bg + num/den


class RamseyStark(object):
    fileData = None
    freqRead = None
    
    timeVec = None
    realVec = None
    imagVec = None
    distVec = None
    
    coeffFit = {"dist": None, "real": None, "imag": None}
    errorFit = {"dist": None, "real": None, "imag": None}
    colorData = {"dist": "green", "real": "red", "imag": "blue"}
    
    def __init__(self, fileData):
        
        self.fileData = fileData
        self.getFreqRead_FromFile()
        self.loadData()
        
        
    def getFreqRead_FromFile(self):
        nameFile = self.fileData.split('\\')[-1]
        strFreq = nameFile.split('_')[-1]
        
        self.freqRead = float(strFreq.partition('G')[0])
        
        
    def getFreqReadGHz(self):
        return self.freqRead
    
    
    def getT2_qubitShift(self):
        mini = +np.inf
        keyBest = None
        for key, error in self.errorFit.items():
            if error[4] < mini:
                mini = error[4]
                keyBest = key
        
        if keyBest is None:
            return [+np.inf, +np.inf, +np.inf, +np.inf]
        
        T2 = self.coeffFit[keyBest][4]
        T2error = self.errorFit[keyBest][4]
        qubitShift = self.coeffFit[keyBest][2]*1e3 / (2*np.pi)
        qubitShifterror = self.errorFit[keyBest][2]*1e3 / (2*np.pi)
        
        return [T2, T2error, qubitShift, qubitShifterror]
        
        
    def loadData(self):
        data = np.loadtxt(self.fileData)
        
        self.timeVec = data[:,0]
        self.realVec = data[:,3]*1000
        self.imagVec = data[:,4]*1000
        self.distVec = data[:,5]*1000
        
        
    def fitSignal(self, key):
        signal = getattr(self, "{}Vec".format(key))
        
        bg = (np.max(signal) + np.min(signal)) / 2.0
        osciAmp = (np.max(signal) - np.min(signal)) / 2.0
        phi0, decayTime = 0, 1000.0
        
        # Wr = 2*np.pi*nbOsci(signal, bg) / self.timeVec[-1]
        Wr = 2*np.pi / Periode(self.timeVec, signal, bg)
                
        guess = (bg, osciAmp, Wr, phi0, decayTime)
        
        bounds_min = (-np.inf, 0, 0, 0, 0)
        bounds_max = (np.max(signal), np.max(signal) - np.min(signal), +np.inf, 2*np.pi, 10000)
        bounds = (bounds_min, bounds_max)
        
        try:
            popt, pcov = curve_fit(ramseyFit, self.timeVec, signal, p0=guess, bounds=bounds)#, maxfev=5000)
        except:
            return (0, 0, 0, 0, 0), (+np.inf, +np.inf, +np.inf, +np.inf, +np.inf)
                    
        return popt, np.sqrt(np.diag(pcov))
    
        
    def fitData(self):
        for key in self.coeffFit.keys():
            coeff, error = self.fitSignal(key)

            self.coeffFit[key] = coeff
            self.errorFit[key] = error
        
    
    def plotterData(self):
        plt.figure('RamseyStark', figsize=(10,6))
        markerStyleData = dict(linestyle='-', marker='',
                           markersize=5, markerfacecoloralt='red', fillstyle='none')
        
        spaceTime = np.linspace(self.timeVec[0], self.timeVec[-1], 5000)
        
        title = ""
        for key in ["dist"]:# self.coeffFit.keys():
            coeffFit = self.coeffFit[key]
            signal = getattr(self, "{}Vec".format(key))
            title += "{}: {}\n".format(key, self.coeffFit[key])
            
            plt.plot(self.timeVec, signal, color=self.colorData[key], **markerStyleData)
            plt.plot(spaceTime, ramseyFit(spaceTime, *coeffFit), ls='-', color=self.colorData[key], lw=2)

        plt.xlabel(r'Waiting time $\tau$ $[\mu s]$')
        plt.ylabel(r'Distance $D_{e-g}$ [mV]')
        plt.title(title)
        
        plt.grid(b=True, color='grey', linestyle=':', linewidth=0.8)
        plt.tight_layout()
        
        plt.show()

        

class NbPhoton(object):
    # meta data
    fileData = None
    rudat = None
    fitted = False
    dictBounds = None
    
    # main data
    arrayFreq = None
    arrayT2 = None
    arrayT2_error = None
    arrayQbShift = None
    arrayQbShift_error = None
    arrayGamma = None
    arrayGamma_error = None
    
    #expected parameters in radian
    dictParam_rad = {
        "bg_qbShift": None, # additive constant for qubit shift
        "bg_gamma": None,   # additive constant for gamma (dephasing)
        "wc": None,         # center frequency
        "k" : None,         # width of resonator
        "chi": None,        # shift of resonator freq by qubit state
        "ed": None         # power on input of cavity ed == |ed|**2
    }

    #fit parameters
    fitParamGamma_rad = None
    fitParamQbShift_rad = None

	# extracted n-phot
    n_phot = []
    n_phot_fit = []
    n_of_phot = None		# We think this is a measured Number of Photons (maximum of n_phot_fit)

    def __init__(self, fileData, rudat, dictBounds, typeLoad, **kwargs_MHz):
        self.fileData = fileData
        self.rudat = rudat
        self.dictBounds = dictBounds
        
        for key in self.dictParam_rad.keys():
            self.dictParam_rad[key] = MHzToRad(kwargs_MHz[key])
        self.dictParam_rad["bg_gamma"] = kwargs_MHz["bg_gamma"]
        
        if typeLoad == "raw":
            self.loadRawData()
        elif typeLoad == "processed":
            self.loadProcessedData()
        else:
            self.loadRawData()
      
    def checkT2_QubitShift(self, T2, qubitShift):
        min_qshift, max_qshift = self.dictBounds["qubitShift"]
        min_t2, max_t2 = self.dictBounds["T2"]
        min_gamma, max_gamma = self.dictBounds["gamma"]
        
        gamma = 1000.0/T2
        if (qubitShift> min_qshift) and (qubitShift< max_qshift):
            if(T2>min_t2) and (T2<max_t2):
                if(gamma>min_gamma) and (gamma<max_gamma):
                    return True
        return False
    
    
    def checkT2_QubitShift_error(self, T2, T2error, qubitShift, qubitShifterror, errorMax = 1):
        
        if T2error < errorMax*T2 and qubitShifterror < errorMax*qubitShift:
            return True
        
        return False
        
        
    def loadRawData(self):
        listFolder = os.listdir(self.fileData)
        
        listFolderRawData = []
        for folder in listFolder:
            if "Ramsey_stark" in folder:
                if "OSC_fit" not in folder:
                    listFolderRawData.append(folder + '\\')
        
        listFileData = []
        for folder in listFolderRawData:
            listFile = os.listdir(self.fileData + folder)
            for file in listFile:
                if file.endswith(".dat"):
                    listFileData.append(self.fileData + folder + file)
            
            
        self.arrayFreq = []
        self.arrayT2 = []
        self.arrayT2_error = []
        self.arrayQbShift = []
        self.arrayQbShift_error = []
        
        nbLoad = len(listFileData)
        sys.stdout.write('Loading data rudat {}dB: {:03d}/{}'.format(self.rudat, 0, nbLoad))
        for i, fileData in enumerate(listFileData):   
            RS = RamseyStark(fileData)
            RS.fitData()
            freqRead = RS.getFreqReadGHz()*1000
            T2, T2error, qubitShift, qubitShifterror = RS.getT2_qubitShift()
            if self.checkT2_QubitShift(T2, qubitShift) and self.checkT2_QubitShift_error(T2, T2error, qubitShift, qubitShifterror):
                self.arrayFreq.append(freqRead)
                self.arrayT2.append(T2)
                self.arrayT2_error.append(T2error)
                self.arrayQbShift.append(qubitShift)
                self.arrayQbShift_error.append(qubitShifterror)
            sys.stdout.write('\rLoading data rudat {}dB: {} / {}'.format(self.rudat, i+1, nbLoad))
            
        sys.stdout.write('\n')
        self.arrayFreq = np.array(self.arrayFreq)
        self.arrayT2 = np.array(self.arrayT2)
        self.arrayT2_error = np.array(self.arrayT2_error)
        self.arrayQbShift = np.array(self.arrayQbShift)
        self.arrayQbShift_error = np.array(self.arrayQbShift_error)
        self.arrayGamma = 1000.0/self.arrayT2 #1000.0 to be in [Mhz] (T2 in [ns])
        self.arrayGamma_error = 1000.0*self.arrayT2_error/self.arrayT2**2
        
            
    def loadProcessedData(self):
        listFolder = os.listdir(self.fileData)
        for folder in listFolder:
            if "Qubit_Shift_vs_Cav_Detun" in folder:
                folderData = folder + '\\'
                break
        
        listFile = os.listdir(self.fileData + folderData)
        for file in listFile:
            if file.endswith(".dat"):
                fileData = file
        
        data = np.loadtxt(self.fileData + folderData + fileData)

              
        self.arrayFreq = []
        self.arrayT2 = []
        self.arrayT2_error = []
        self.arrayQbShift = []
        self.arrayQbShift_error = []
        
        for i, freqRead in enumerate(data[:,0]):   
            T2, T2error, qubitShift, qubitShifterror = data[i,1], data[i,2], data[i,3], data[i,4]
            if self.checkT2_QubitShift(T2, qubitShift) and self.checkT2_QubitShift_error(T2, T2error, qubitShift, qubitShifterror):
                self.arrayFreq.append(freqRead*1000)
                self.arrayT2.append(T2)
                self.arrayT2_error.append(T2error)
                self.arrayQbShift.append(qubitShift)
                self.arrayQbShift_error.append(qubitShifterror)
            
        self.arrayFreq = np.array(self.arrayFreq)
        self.arrayT2 = np.array(self.arrayT2)
        self.arrayT2_error = np.array(self.arrayT2_error)
        self.arrayQbShift = np.array(self.arrayQbShift)
        self.arrayQbShift_error = np.array(self.arrayQbShift_error)
        self.arrayGamma = 1000.0/self.arrayT2 #1000.0 to be in [Mhz] (T2 in [ns])
        self.arrayGamma_error = 1000.0*self.arrayT2_error/self.arrayT2**2

        
        
    def getGuess(self, typeBg):
        listKey = ["bg_{}".format(typeBg), "wc", "k", "chi", "ed"]
        limMini = [-np.inf, 0, 0, 0, 0]
        limMax = [+np.inf, 2, 2, 2, +np.inf]
        
        guess = [self.dictParam_rad[key] for key in listKey]
        bounds_min = ([mini*x for mini, x in zip(limMini, guess)])
        bounds_max = ([maxi*x for maxi, x in zip(limMax, guess)])
        
        return guess, (bounds_min, bounds_max)
        
        
    def fitGamma(self):
        
        guess, bounds = self.getGuess("gamma")
        xData = MHzToRad(self.arrayFreq)
        yData, yError = self.arrayGamma, self.arrayGamma_error
        
        try:
            popt, _ = curve_fit(ramseyStarkDephasing, xData, yData, sigma=yError, p0=guess, bounds=bounds)
            self.fitParamGamma_rad = popt
        except:
            popt, _ = curve_fit(ramseyStarkDephasing, xData, yData, sigma=yError, p0=guess, maxfev=5000)
            self.fitParamGamma_rad = popt
        
        
    def fitQbShift(self):
        
        guess, bounds = self.getGuess("qbShift")
        xData = MHzToRad(self.arrayFreq)
        yData, yError = MHzToRad(self.arrayQbShift), MHzToRad(self.arrayQbShift_error)
        try:
            popt, _ = curve_fit(ramseyStarkQbShift, xData, yData, sigma=yError, p0=guess, bounds=bounds)
            self.fitParamQbShift_rad = popt
        except:
            popt, _ = curve_fit(ramseyStarkQbShift, xData, yData, sigma=yError, p0=guess, maxfev=5000)
            self.fitParamQbShift_rad = popt
                
        
    def fitGamma_QbShift(self):
        guessGamma = self.getGuess("gamma")
        guessQbShift = self.getGuess("qbShift")
        guess = [guessGamma[0]] + guessQbShift
        
        xData = MHzToRad(self.arrayFreq)
        yGamma, yGammaError = self.arrayGamma, self.arrayGamma_error
        yQbShift, yQbShiftError = MHzToRad(self.arrayQbShift), MHzToRad(self.arrayQbShift_error)
        
        comboX = np.append(xData, xData)
        comboY = np.append(yGamma, yQbShift)
        comboYerror = np.append(yGammaError, yQbShiftError)
        
        def comboFunc(comboData, n, bg_gamma, bg_qbShift, *param):
            # single data set passed in, extract separate data
            xData1 = comboData[:n//2] # first data
            xData2 = comboData[n//2:] # second data
        
            result1 = ramseyStarkDephasing(xData1, bg_gamma, *param)
            result2 = ramseyStarkQbShift(xData2, bg_qbShift, *param)
        
            return np.append(result1, result2)
        
        func = lambda comboData, bg_gamma, bg_qbShift, *param: comboFunc(comboData, len(xData), bg_gamma, bg_qbShift, *param)
        popt, _ = curve_fit(func, comboX, comboY, sigma=comboYerror, p0=guess)
        self.fitParamQbShift_rad = np.array([popt[1], *popt[2:]])
        self.fitParamGamma_rad = np.array([popt[0], *popt[2:]])
        
        
        
    def fitNbPhoton(self):
        k_f = self.fitParamGamma_rad[2]
        chi_f = self.fitParamGamma_rad[3]
        
        minFreq = np.min(self.arrayFreq)
        maxFreq = np.max(self.arrayFreq)
        freqSpace = np.linspace(minFreq, maxFreq, 5000)
        fitGamma_rad = ramseyStarkDephasing(MHzToRad(freqSpace), *self.fitParamGamma_rad)
        
        self.n_phot = self.arrayGamma *( (k_f**2  +  chi_f**2)/(2 * k_f * chi_f**2) )
        self.n_phot_fit = fitGamma_rad *( (k_f**2  +  chi_f**2)/(2 * k_f * chi_f**2) )

        self.n_of_phot = np.max(self.n_phot_fit)
        
        
    def fitData(self, typeFit):
        if typeFit == "independant":
            self.fitGamma()
            self.fitQbShift()
        else:
            self.fitGamma_QbShift()
        
        self.fitNbPhoton()
        
        self.fitted = True
        
        
    def plotterQbShift(self):
        plt.figure('Qubit shift from Ramsey Stark experiment, rudat {}dB'.format(self.rudat))
        plt.grid()
        
        plt.plot(self.arrayFreq, self.arrayQbShift, '.', label="experimental data")
        plt.errorbar(self.arrayFreq, self.arrayQbShift, yerr=self.arrayQbShift_error, ls='none', color='grey')
        plt.title("Raw data, rudat {}dB".format(self.rudat))
        
        if self.fitted:
            minFreq = np.min(self.arrayFreq)
            maxFreq = np.max(self.arrayFreq)
            freqSpace = np.linspace(minFreq, maxFreq, 5000)
            fitQbShift_rad = ramseyStarkQbShift(MHzToRad(freqSpace), *self.fitParamQbShift_rad)
            plt.plot(freqSpace, RadToMHz(fitQbShift_rad), label='fit')
            labelParam = "Rudat: {}dB, bg: {:.2f}MHz, wc: {:.1f}MHz,\n k: {:.1f}MHz, chi: {:.1f}MHz, ed: {:.1f}MHz".format(self.rudat, *RadToMHz(self.fitParamQbShift_rad))
            plt.title(labelParam)
        
        plt.xlabel('Cavity frequency  [MHz]')
        plt.ylabel('Qubit shift')
        plt.legend()
        plt.show()
        
    def plotterGamma(self):
        plt.figure('1/T2 from Ramsey Stark experiment, rudat {}dB'.format(self.rudat))
        plt.grid()
        
        
        plt.plot(self.arrayFreq, self.arrayGamma, '.', label="experimental data")
        plt.errorbar(self.arrayFreq, self.arrayGamma, yerr=self.arrayGamma_error, ls='none', color='grey')
        plt.title("Raw data, rudat {}dB".format(self.rudat))
        
        if self.fitted:
            minFreq = np.min(self.arrayFreq)
            maxFreq = np.max(self.arrayFreq)
            freqSpace = np.linspace(minFreq, maxFreq, 5000)
            fitGamma_rad = ramseyStarkDephasing(MHzToRad(freqSpace), *self.fitParamGamma_rad) 
            plt.plot(freqSpace, fitGamma_rad, label='fit')
            labelParam = "Rudat: {}dB, bg: {:.2f}MHz, wc: {:.1f}MHz,\n k: {:.1f}MHz, chi: {:.1f}MHz, ed: {:.1f}MHz".format(self.rudat, *RadToMHz(self.fitParamGamma_rad))
            plt.title(labelParam)
        
        plt.xlabel('Cavity frequency  [MHz]')
        plt.ylabel(r'$\Gamma$ [rad]')
        plt.legend()
        plt.show()
        
    def plotterT2(self):
        plt.figure('T2 from Ramsey Stark experiment, rudat {}dB'.format(self.rudat))
        plt.grid()
        
        plt.plot(self.arrayFreq, self.arrayT2, '.', label="experimental data")
        plt.errorbar(self.arrayFreq, self.arrayT2, yerr=self.arrayT2_error, ls='none', color='grey')
        plt.title("Raw data, rudat {}dB".format(self.rudat))
        
        plt.xlabel('Cavity frequency  [MHz]')
        plt.ylabel(r'T2 [ns]')
        plt.legend()
        plt.show()
        
        
    def plotterNbPhoton(self):
        plt.figure('Number of photons, rudat {}dB'.format(self.rudat))
        plt.grid()
        
        plt.plot(self.arrayFreq, self.n_phot, '.', label="experimental data")
        plt.title("Raw data, rudat {}dB".format(self.rudat))
        
        if self.fitted:
            minFreq = np.min(self.arrayFreq)
            maxFreq = np.max(self.arrayFreq)
            freqSpace = np.linspace(minFreq, maxFreq, 5000)
            plt.plot(freqSpace, self.n_phot_fit, '--', label='fit')
            plt.title('Photon number is {:.3f} for rudat {}dB'.format(self.n_of_phot, self.rudat))
            
        plt.xlabel('Cavity frequency  [MHz]')
        plt.ylabel('Number of photons in polariton')
        plt.legend()
        plt.show()

    
def plotterGlobal(listFolder, listRudat, dictBounds, typeLoad, typeFit, **kwargs_MHz):

    cmap = plt.get_cmap('jet')
    listColor = [cmap(i) for i in np.linspace(0, 1.0, len(listFolder))]
    
    figPhoton = plt.figure('Number of photons', figsize=[10,8])
    axPhoton = figPhoton.add_subplot(111)    
    axPhoton.grid()
    axPhoton.set_xlabel('Cavity frequency  [MHz]')
    axPhoton.set_ylabel('Number of photons in cavity')
    
    figGamma = plt.figure('1/T2 from Ramsey Stark experiment', figsize=[10,8])
    axGamma = figGamma.add_subplot(111)    
    axGamma.grid()
    axGamma.set_xlabel('Cavity frequency  [MHz]')
    axGamma.set_ylabel(r'$\Gamma$ [rad]')
    
    figQbShift = plt.figure('Qubit shift from Ramsey Stark experiment', figsize=[10,8])
    axQbShift = figQbShift.add_subplot(111)    
    axQbShift.grid()
    axQbShift.set_xlabel('Cavity frequency  [MHz]')
    axQbShift.set_ylabel('Qubit shift')
    
    figMaxPhoton = plt.figure('Cavity photon number', figsize=[10,8])
    axMaxPhoton = figMaxPhoton.add_subplot(111)    
    axMaxPhoton.grid()
    axMaxPhoton.set_xlabel('Rudat attenuation [dB]')
    axMaxPhoton.set_ylabel('Maximum number of photons in cavity')
    
    listNbPhoton = []
    for folder, rudat, color in zip(listFolder, listRudat, listColor):
        ACs = NbPhoton(folder, rudat, dictBounds[rudat], typeLoad=typeLoad, **kwargs_MHz)
        ACs.fitData(typeFit)
        listNbPhoton.append(ACs.n_of_phot)
        
        minFreq = np.min(ACs.arrayFreq)
        maxFreq = np.max(ACs.arrayFreq)
        freqSpace = np.linspace(minFreq, maxFreq, 5000)
        fitGamma_rad = ramseyStarkDephasing(MHzToRad(freqSpace), *ACs.fitParamGamma_rad) 
        fitQbShift_rad = ramseyStarkQbShift(MHzToRad(freqSpace), *ACs.fitParamQbShift_rad)
        
        axPhoton.plot(ACs.arrayFreq, ACs.n_phot, '.', color=color, label="rudat {}dB - n={:.3f}".format(ACs.rudat, ACs.n_of_phot))
        axPhoton.plot(freqSpace, ACs.n_phot_fit, color=color)
        
        axGamma.plot(ACs.arrayFreq, ACs.arrayGamma, '.', color=color, label="rudat {}dB".format(ACs.rudat))
        # axGamma.errorbar(ACs.arrayFreq, ACs.arrayGamma, yerr=ACs.arrayGamma_error, ls='none', color='grey')
        axGamma.plot(freqSpace, fitGamma_rad, color=color)
    
        axQbShift.plot(ACs.arrayFreq, ACs.arrayQbShift, '.', color=color, label="rudat {}dB".format(ACs.rudat))
        # axQbShift.errorbar(ACs.arrayFreq, ACs.arrayQbShift, yerr=ACs.arrayQbShift_error, ls='none', color='grey')
        axQbShift.plot(freqSpace, RadToMHz(fitQbShift_rad), color=color)
        
        axMaxPhoton.plot( [rudat], [ACs.n_of_phot], '.', markersize=15, color=color, label="rudat {}dB - n={:.3f}".format(ACs.rudat, ACs.n_of_phot))
        
        
    funcFit = lambda rudat, a, b: a*rudat + b
    popt, _ = curve_fit(funcFit, listRudat, listNbPhoton, p0=(-0.01, 0.5))
    print("Coefficient for the linear fir", popt)
    axMaxPhoton.plot( listRudat, funcFit(np.array(listRudat), *popt), '-k')
    
    axGamma.legend()
    axQbShift.legend()
    axPhoton.legend()
    axMaxPhoton.legend()
    
    
    
   
##--data----####   
baseFolder = "D:\\data_12_01_19\\Documents\\Stage_Neel\\All_data\\Nb_photon\\After_confinement\\"
    

R27_17 = baseFolder + "17062020_Rudat_27dB\\"
R28_11 = baseFolder + "11062020_Rudat_28dB\\"
R28_17 = baseFolder + "17062020_Rudat_28dB\\"
R29_17 = baseFolder + "17062020_Rudat_29dB\\"
R30_10 = baseFolder + "10062020_Rudat_30dB\\"
# R30_17 = baseFolder + "17062020_Rudat_30dB\\"


# can be "raw" or "processed"
#if raw : redo the fit of raw data
#if processed : used value from the fit already made
typeLoad = "raw" 

# can be "independant" or "together"
#if independant : make the fit of Gamma and Qubit shift independantly
#if together : make one fit for Gamma and Qubit shift together
typeFit = "independant" 

##--expectations----####
freqQb 	   = 6.29911            #[GHz] #value of qubit freq, calibrated with Ramsey
freq_used  = 6.29911+0.01		#[GHz] #value of setted freq_qubit
deltaFreq  = freqQb - freq_used
T2_ns	   = 4500.0             # average T2 measured far from resonance without pump


kwargs_MHz = {
    "bg_qbShift":   abs(deltaFreq)*1e-3,  # additive constant for qubit shift
    "bg_gamma":     1000.0/T2_ns,    # additive constant for gamma (dephasing)
    "wc":           7074.5,          # center frequency
    "k" :           10.00,           # width of resonator
    "chi":          7.75,            # shift of resonator freq by qubit state
    "ed":           10              # power on input of cavity ed == |ed|**2
    }
    

#min and max value of different quantities
dictBoundsData = {27: {"qubitShift": (9.3, 10.2), "T2": (10, 20000), "gamma": (0,10)},
                  28: {"qubitShift": (9.3, 10.2), "T2": (10, 20000), "gamma": (0,8)},
                  29: {"qubitShift": (9.3, 10.2), "T2": (10, 20000), "gamma": (0,6.5)},
                  30: {"qubitShift": (9.3, 10.2), "T2": (10, 20000), "gamma": (0,5.2)}}


plt.close('all')


######################################
####### Single file processing #######
######################################

rudat=27
ramseyStark = NbPhoton(R27_17, rudat, dictBoundsData[rudat], typeLoad=typeLoad, **kwargs_MHz)
ramseyStark.fitData(typeFit)
ramseyStark.plotterQbShift()
ramseyStark.plotterGamma()
if ramseyStark.fitted:
    ramseyStark.plotterNbPhoton()
else:
    ramseyStark.plotterT2()
   
    
   
########################################
####### Multiple file processing #######
########################################

# listFolder = [R30_10, R29_17, R28_17, R27_17]
# listRudat = [30,29,28,27]
# plotterGlobal(listFolder, listRudat, dictBoundsData, typeLoad=typeLoad, typeFit=typeFit, **kwargs_MHz)

   
plt.show()
