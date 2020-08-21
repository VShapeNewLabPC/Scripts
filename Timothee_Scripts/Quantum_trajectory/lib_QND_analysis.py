import os
import gc
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

MEDIUM_SIZE = 19
BIGGER_SIZE = 21

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


dark = False
if dark:
    import matplotlib as mpl
    mpl.rc('xtick', color='w')
    mpl.rc('ytick', color='w')
    mpl.rc('grid', color='w')


def fitSmoothingRegime(t, a, b):
    return a*np.exp( -b*t )

def fitSmoothingRegimeLabel(var, a, b):
    if b>0:
        return r'$%.3fexp[-%.2e\times %s]$'%(a, b, var)
    elif b==0:
        return r'$%.3f$'%(a)
    else:
        return r'$%.3fexp[%.2e\times %s]$'%(a, abs(b), var)


############################################################################
###################   CLASS QNDanalyse   ###################################
############################################################################  

class QNDanalyse(object):
    nameQNDfolder = None
    with_pio2 = None
    
    listAllNaver = None
    dictNaver = None
    dictLenNaver = None
    
    listRudat = None
    lenRudat = None
    
    dictFileTraject = None
    dictData = None
    dictCoeffFitJR = None
    
    listKeyQND = ['Q', 'Qe', 'Qg', 'uncertainty', 
                    'Pgg', 'Pee', 'Peg', 'Pge']
    listKeyJump = ['n_jump', 'n_not_jump', 
                     'decay_rate_', 'excitation_rate_',
                     'decay_rate_censor', 'excitation_rate_censor',
                     'mean_duration_e', 'mean_duration_g', 
                     'mean_duration_e_censor', 'mean_duration_g_censor', 
                     'decay_rate_proba', 'excitation_rate_proba', 
                     'decay_rate_stat', 'excitation_rate_stat']
    listKeyData = listKeyQND + listKeyJump
    
    listKeyPlot = None
    listColor = None
    
    def __init__(self, list_QND_folder, listRudat, dictNaver, with_pio2, **kwargs_compute):
        start = time.time()
        sys.stdout.write('Initialisation of object ...')
        
        self.with_pio2 = with_pio2
            
        self.setNameQNDFolder(**kwargs_compute)
        self.setKeyPlot(**kwargs_compute)
        
        self._init_Rudat(listRudat)
        self._init_Naver(dictNaver, list_QND_folder)

        self.dictCoeffFitJR = {}
        for rudat in self.listRudat:
            self.dictCoeffFitJR[rudat] = {}
        
        self.dictFileTraject = {}
        name_folder = 'QNDness_' + self.nameQNDfolder + '\\'
        for QND_folder, rudat in zip(list_QND_folder, self.listRudat):
            
            self.dictFileTraject[rudat] = {}                
            for Naver in self.dictNaver[rudat]:
                if type(rudat) is str:
                    name_file = 'QNDness_R{}_naver{}.dat'.format(0, Naver)
                else:
                    name_file = 'QNDness_R{}_naver{}.dat'.format(rudat, Naver)
                self.dictFileTraject[rudat][Naver] = QND_folder + name_folder + name_file

        self.dictData = {}
        for index, key in enumerate(self.listKeyData):
            self.dictData[key] = {}
            for rudat in self.listRudat:
                self.dictData[key][rudat] = np.zeros((self.dictLenNaver[rudat]))
            
        end = time.time()
        sys.stdout.write('\rInitialisation of object ... done - {:.3f}s'.format(end-start))


    def _init_Rudat(self, listRudat):
        self.listRudat = listRudat
        self.lenRudat = len(self.listRudat)


    def _init_Naver(self, dictNaver, list_QND_folder):
        self.dictNaver = {}
        self.dictLenNaver = {}
        self.listAllNaver = []
        
        name_folder = 'QNDness_' + self.nameQNDfolder + '\\'
        for QND_folder, rudat in zip(list_QND_folder, self.listRudat):
            
            path = QND_folder + name_folder
            listdir = os.listdir(path)
            listAllNaver = sorted([self.getNaver_fromFile(file) for file in listdir])
            
            try:
               allFlag = dictNaver[rudat]['all'] 
            except KeyError:
                allFlag = True
                
            if allFlag:
                self.listAllNaver = self.listAllNaver + listAllNaver
                self.dictNaver[rudat] = np.array(listAllNaver)
                self.dictLenNaver[rudat] = len(listAllNaver)
            else:
                startNaver = dictNaver[rudat]['start']
                if startNaver is None:
                    sartNaver = np.min(listAllNaver)
                    
                stopNaver = dictNaver[rudat]['stop']
                if stopNaver is None:
                    stopNaver = np.max(listAllNaver)
                    
                exceptNaver = dictNaver[rudat]['except']
                
                listNaver = []
                for Naver in listAllNaver:
                    if Naver in exceptNaver:
                        continue
                    
                    if startNaver <= Naver and Naver <=stopNaver:
                        self.listAllNaver.append(Naver)
                        listNaver.append(Naver)
                  
                
                self.dictNaver[rudat] = np.array(listNaver)
                self.dictLenNaver[rudat] = len(listNaver)
                
        self.listAllNaver = np.unique(self.listAllNaver)
            
    
    def loadQNDresult(self):
        start = time.time()
        
        print('')
        for rudat in self.listRudat:
            for ind_naver, Naver in enumerate(self.dictNaver[rudat]):
                sys.stdout.write('\rLoad QND result of rudat={}dB and naver={}ns ...'.format(rudat, Naver))
                
                file = self.dictFileTraject[rudat][Naver]
                try:
                    data = np.loadtxt(file)
                except OSError:
                    oldNameFile = 'QNDness_R{}_naver{}.dat'.format(rudat, Naver)
                    newNameFile = 'QNDness_R0_naver{}.dat'.format(Naver)
                    file = file.replace(oldNameFile, newNameFile)
                    data = np.loadtxt(file)
                
                for index, key in enumerate(self.listKeyData):
                    self.dictData[key][rudat][ind_naver] = data[index]
           
        
        for keyType, keyState in [('decay', 'e'), ('excitation', 'g')]:
            for keyData in ['', '_censor']:
                key = '{}_rate_mean{}'.format(keyType, keyData)
                self.dictData[key] = {}
                for rudat in self.listRudat:
                    keyOld = 'mean_duration_{}{}'.format(keyState, keyData)
                    dataOld = self.dictData[keyOld][rudat]
                    
                    self.dictData[key][rudat] = 1/dataOld

        end = time.time()
        sys.stdout.write('\rLoading QND result ... done - {:.3f}s                             \n'.format(end-start))
       
        
    def setNameQNDFolder(self, **kwargs):
        definition = kwargs['definition'] #'Me' or 'Sl'
        threshold = kwargs['threshold'] #'mean' or 'fit' or 'double'
        security = kwargs['security'] #True or False
        
        self.nameQNDfolder = ''
        if threshold == 'double':
            self.nameQNDfolder += 'double_{}'.format(definition)
        else:
            self.nameQNDfolder += 'single_{}'.format(threshold)
            
        if not security:
            self.nameQNDfolder += '_bypass'
      
        
    def setKeyPlot(self, **kwargs_compute):
        dataToPlot = kwargs_compute['dataToPlot']
        
        self.listColor = []
        self.listKeyPlot = []
        if 'rate' in dataToPlot:
            self.listColor += ['tab:purple', 'orange']
            self.listKeyPlot += ['decay_rate_', 'excitation_rate_']
            
        if 'rate_censor' in dataToPlot:
            self.listColor += ['tab:blue', 'tab:red']#['tab:red']#['tab:blue', 'tab:red']
            self.listKeyPlot += ['decay_rate_censor', 'excitation_rate_censor']#['excitation_rate_censor']#['decay_rate_censor', 'excitation_rate_censor']
            
        if 'rate_proba' in dataToPlot:
            self.listColor += ['darkviolet', 'darkorange']
            self.listKeyPlot += ['decay_rate_proba', 'excitation_rate_proba']
            
        if 'rate_stat' in dataToPlot:
            self.listColor += ['darkblue', 'darkred']
            self.listKeyPlot += ['decay_rate_stat', 'excitation_rate_stat']
            
            
    def getLabelRate(self, typeRate):
        listType = typeRate.split('_')
        return ' '.join(listType)
            
    
    def getSavePath(self, list_QND_folder, QND_folder_base):
        
        if self.lenRudat == 1 and self.listRudat[0] == 0:
            savepath = list_QND_folder[0] + 'QND_vs_Rudat\\' + self.nameQNDfolder + '\\' 
        else:
            savepath = QND_folder_base + 'QND_vs_Rudat\\' + self.nameQNDfolder + '\\'
        
        return savepath
        
        
    ##############################
    ###   GET FUNCTIONS  #########
    def getRudat_fromFolder(self, folder):
        ind = folder.find('R')
        if ind == -1:
            return 0
        
        str_rudat = folder[ind+1:ind+3]
        return int(str_rudat)        
    
    def getNaver_fromFile(self, file):
        end = file.partition('naver')[2]
        str_Naver = end.partition('.')[0]
        
        return int(str_Naver)
    
    def getListDataVsRudat(self, key, Naver):
        listRudat = []
        listData = []
        for rudat in self.listRudat:
            if Naver in self.dictNaver[rudat]:
                indNaver = np.argwhere(self.dictNaver[rudat]==Naver)[0,0]
                
                listRudat.append(rudat)
                listData.append(self.dictData[key][rudat][indNaver])
                
        return np.array(listData), np.array(listRudat)
    
    
    def getMinMax_Q(self, Naver):
        listQ, _ = self.getListDataVsRudat('Q', Naver)
        listQe, _ = self.getListDataVsRudat('Qe', Naver)
        listQg, _ = self.getListDataVsRudat('Qg', Naver)
        
        return (np.min([np.min(listQ), np.min(listQe), np.min(listQg)]),
                np.max([np.max(listQ), np.max(listQe), np.max(listQg)]))
    
    
    def setListNaver_fromMaxQe(self):
        listMaxi = []
        self.maxQeNaver = []
        for rudat in self.listRudat:
            indmax = np.argmax(self.dictData['Q'][rudat])
            maxi = np.max(self.dictData['Q'][rudat])
            Naver = self.dictNaver[rudat][indmax]
            
            listMaxi.append(maxi)
            self.maxQeNaver.append(Naver)
            
            print("Rudat={}dB, max QND-ness={:.1f}, Naver={}".format(rudat, maxi*100, Naver))


    
    #######################################
    ###   FIT Jump Rate VS Naver  #########    
    def getNaverBoundary(self, rudat, dtype, dictNaverFit):
        listType = dtype.split('_')
        gammaType = listType[0]
        try:
            startNaver = dictNaverFit[gammaType][rudat]['start']
            stopNaver  = dictNaverFit[gammaType][rudat]['stop']
            uncertainty = dictNaverFit[gammaType][rudat]['uncertainty']
        except KeyError:
            startNaver = None
            stopNaver  = None
            uncertainty = None
        
        return startNaver, stopNaver, uncertainty
    
    
    def getIndBoundary(self, rudat, startNaver, stopNaver):
        if startNaver is None:
            startNaver = self.dictNaver[rudat][0]
        if stopNaver is None:
            stopNaver = self.dictNaver[rudat][-1]
            
        return (np.argwhere(self.dictNaver[rudat] == startNaver)[0, 0],
                np.argwhere(self.dictNaver[rudat] == stopNaver)[0, 0])
      
        
    def fitJumpRate_SmoothingRegime(self, rudat, dtype, dictNaverFit):    
        startNaver, stopNaver, uncertainty = self.getNaverBoundary(rudat, dtype, dictNaverFit)
        indStart, indStop = self.getIndBoundary(rudat, startNaver, stopNaver)
        
        listNaver = self.dictNaver[rudat][indStart:indStop]
        listData = self.dictData[dtype][rudat][indStart:indStop]*1000
                    
        sigma = np.ones(np.shape(listNaver))
        if uncertainty is not None:
            indStart, indStop = self.getIndBoundary(rudat, uncertainty[0], uncertainty[1])
            sigma[indStart:indStop] = uncertainty[2]
            
        a0, a_min, a_max = 1, 0, 100
        b0, b_min, b_max = 0, -np.inf, np.inf
        
        p0 = (a0, b0)
        bounds_min = (a_min, b_min)
        bounds_max = (a_max, b_max)
    
        coeffFit, error = curve_fit(fitSmoothingRegime, listNaver, listData, 
                                     p0=p0, bounds=(bounds_min, bounds_max),
                                     sigma=sigma)
    
        self.dictCoeffFitJR[rudat][dtype] = coeffFit
        
    
    def fitAllJumpRate(self, dictNaverFit):
        
        start = time.time()
        sys.stdout.write('Fit of jump rates ...')
                
        for rudat in self.listRudat:
            
            for dtype in self.listKeyPlot:
                self.fitJumpRate_SmoothingRegime(rudat, dtype, dictNaverFit)
                
        end = time.time()
        sys.stdout.write('\rFit of jump rates ...   done - {:.3f}s\n'.format(end-start))

            
    #############################
    ###   PLOT QNDness  #########
    def plotQNDnessVsRudat(self, Naver, savepath='', savename='QNDness_vs_rudat.png'):        
        listQ, listRudat   = self.getListDataVsRudat('Q', Naver)
        listUncertainty, _ = self.getListDataVsRudat('uncertainty', Naver)
        listQe, _ = self.getListDataVsRudat('Qe', Naver)
        listQg, _ = self.getListDataVsRudat('Qg', Naver)
        
        fig = plt.figure(num=0, figsize=[18,8])
        fig.tight_layout()
        
        list_mini, list_maxi = [], []
        for N in self.listAllNaver:
            mini, maxi = self.getMinMax_Q(N)
            list_mini.append(mini)
            list_maxi.append(maxi)
        ymin, ymax = np.min(list_mini), np.max(list_maxi)      
        
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
        ax.set_ylim(ymin-0.05, ymax+0.05)
        
        errorBar = np.array(listUncertainty) / 100
        ax.errorbar(listRudat, listQ, yerr=errorBar, fmt='o-g', linewidth=1.5, label='Q')
        ax.plot(listRudat, listQe, 'o-r', linewidth=1.5, label='Qe')
        ax.plot(listRudat, listQg, 'o-b', linewidth=1.5, label='Qg')
        
        ax.set_title('QNDness vs rudat for $|e\\rangle$ and $|g\\rangle$ state for a time averaging of {}ns'.format(Naver))
        ax.set_xlabel('Rudat (dB)')
        ax.set_ylabel('QNDness')
            

        fig.legend()
        if savepath != '':
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            plt.savefig(savepath+savename)
            plt.cla()
            plt.clf()
            plt.close(fig)
            gc.collect()
        else:
            plt.show()


    def plotQNDnessVsNaver(self, rudat, savepath='', savename='QNDness_vs_naver.png'):
        listQ = self.dictData['Q'][rudat]
        listUncertainty = self.dictData['uncertainty'][rudat]/100
        listQe = self.dictData['Qe'][rudat]
        listQg = self.dictData['Qg'][rudat]
        
        
        fig = plt.figure(num=1, figsize=[18,8])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
        
        ax.errorbar(self.dictNaver[rudat], listQ, yerr=listUncertainty, 
                    fmt='o-g', linewidth=1.5, label='Q')
        
        ax.plot(self.dictNaver[rudat], listQe, 'o-r', linewidth=1.5, label='Qe')
        ax.plot(self.dictNaver[rudat], listQg, 'o-b', linewidth=1.5, label='Qg')
        ax.axvline(x=self.dictNaver[rudat][np.argmax(listQe)], color='black')
        
        ax.set_title('QNDness versus time averaging for $|e\\rangle$ and $|g\\rangle$ state with rudat={}dB'.format(rudat))
        ax.set_xlabel('Time averaging [ns]')
        ax.set_ylabel('QNDness')

        ax.legend()
        if savepath != '':
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            plt.savefig(savepath+savename)
            plt.cla()
            plt.clf()
            plt.close(fig)
            gc.collect()
        else:
            plt.show()
        
    
    
    #############################
    ###   PLOT JumpRate  ########
    def plotJumpRateVsRudat(self, Naver, savepath='', savename='Jump_rate_vs_rudat.png'):        
        
        fig = plt.figure(num=2, figsize=[18,8])
        fig.tight_layout()
        
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
        
        listType = self.listKeyPlot
        listColor = self.listColor
        for dtype, color in zip(listType, listColor):

            labelType = self.getLabelRate(dtype)
            listData, listRudat = self.getListDataVsRudat(dtype, Naver)
            ax.plot(listRudat, listData*1000, 'o-', color=color, linewidth=1.5, label=labelType)
            
            
        
        ax.set_title('Relaxation and excitation rate vs rudat for a time averaging of {}ns'.format(Naver))
        ax.set_xlabel('Rudat (dB)')
        ax.set_ylabel('Jump rate [$µs^{-1}$]')
        
        fig.legend()
        if savepath != '':
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            plt.savefig(savepath+savename)
            plt.cla()
            plt.clf()
            plt.close(fig)
            gc.collect()
        else:
            plt.show()
    
    
    def plotJumpRateVsNaver(self, rudat, fit=True, savepath='', savename='Jump_rate_vs_naver'):
        
        fig = plt.figure(num=3, figsize=[18,8])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
    
        listType = self.listKeyPlot
        listColor = self.listColor
        
        listMin = []
        listMax = []

        NaverSpace = np.linspace(self.dictNaver[rudat][0], self.dictNaver[rudat][-1], 5000)
        for dtype, color in zip(listType, listColor):
            
            listData = self.dictData[dtype][rudat]*1000
            listMin.append(np.min(listData))
            listMax.append(np.max(listData))
            
            labelType = self.getLabelRate(dtype)
            ax.plot(self.dictNaver[rudat], listData, 'o', color=color, linewidth=1.5, label=labelType)
            
            if fit:
                coeffFit = self.dictCoeffFitJR[rudat][dtype]
                data_fit = fitSmoothingRegime(NaverSpace, *coeffFit)
                ax.plot(NaverSpace, data_fit, '--', color=color)#, label=label)

                # if dtype in ['decay_rate_censor', 'excitation_rate_censor']:
                #     data_fit = fitSmoothingRegime(NaverSpace, *coeffFit)
                #     # label = labelType + ' fit: ' + fitSmoothingRegimeLabel('N_{aver}', *coeffFit)
                # else:
                #     ax.axhline(y=coeffFit, ls='--', color=color)#, label=label)
                
                

        ax.set_ylim(0.95*np.min(listMin), 1.05*np.max(listMax)) 
        ax.set_yscale('log')
        ax.set_title('Relaxation and excitation rate vs time averaging with rudat={}dB'.format(rudat))
        ax.set_xlabel('Number of averaging')
        ax.set_ylabel('Jump rate [$µs^{-1}$]')

        ax.legend()
        if savepath != '':
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            plt.savefig(savepath+savename)
            plt.cla()
            plt.clf()
            plt.close(fig)
            gc.collect()
        else:
            plt.show()
            
    
    def plotJumpRateFitVsRudat(self, tread, dictNaverFit, RudatToNbPhoton=None, savepath=''):
        fig = plt.figure(figsize=(13,7))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.grid()
            
        listType = self.listKeyPlot
        listColor = self.listColor
        dictJumpRate = {}
        for dtype, color in zip(listType, listColor):
            
            if "excitation" in dtype:
                color = "red"
            else:
                color="blue"
            if "censor" in dtype:
                marker = "P"
                ls = ':'
            else:
                marker = "s"
                ls = '-'
            
            dictJumpRate[dtype] = []
            for rudat in self.listRudat:
                coeffFit = self.dictCoeffFitJR[rudat][dtype]
                dictJumpRate[dtype].append(fitSmoothingRegime(0, *coeffFit))  
            
            ax1.plot(self.listRudat, dictJumpRate[dtype], ls=ls, marker=marker, color=color)
            
        ax1.plot([], [], ls=':', marker="P", color='grey', label='Exponential decaying')
        ax1.plot([], [], ls=':', marker="s", color='grey', label='Master equation')
            
        # ax1.set_ylim(0, 1.1*np.max(list_Naver0)) 
        # ax1.set_title('Extracted decay and excitation rates from fit vs rudat')
        ax1.set_ylabel('Jump rate [$µs^{-1}$]')
        ax1.set_yscale('log')
        ax1.legend(loc='upper right')
        ax1.set_xlabel('Rudat attenuation [dB]')
        
        if RudatToNbPhoton is not None:
            ax12 = ax1.twiny()
            R, Nb = RudatToNbPhoton[0], RudatToNbPhoton[1]
            listX = []
            listRudat= [14,16,18,20,22,24]
            for rudat in listRudat:
                coeff = 10**((R-rudat)/10)
                listX.append(Nb*coeff)
            
            ax12.set_xlabel(r'Photon number, $n_{ph}$')
            newTickLocations = np.array(listRudat)
            ax12.set_xlim(ax1.get_xlim())
            ax12.set_xticks(newTickLocations)
            ax12.set_xticklabels(["{:.1f}".format(nb) for nb in listX])
            
        fig.tight_layout()
        
        
        fig2 = plt.figure(figsize=(10,7))
        ax2 = fig2.add_subplot(1, 1, 1)
        ax2.grid() 
        
        sumJR = np.array(dictJumpRate['decay_rate_censor']) + \
                np.array(dictJumpRate['excitation_rate_censor'])
        ax2.plot(self.listRudat, 1000/sumJR, marker = "P", ls = ':', color='#FF6C00', label='Exponential decaying')
        sumJR = np.array(dictJumpRate['decay_rate_proba']) + \
                np.array(dictJumpRate['excitation_rate_proba'])
        ax2.plot(self.listRudat, 1000/sumJR, marker = "s", ls = '-', color='green', label='Master equation')
        ax2.axhline(y=3220, color="black", lw=3, ls='--', label=r'Relaxation $T_1$')
        ax2.set_xlabel('Rudat attenuation [dB]')
        ax2.set_ylabel('Relaxation time $T_1$ [ns]')
        ax2.legend(loc='upper right')
        
        if RudatToNbPhoton is not None:
            ax22 = ax2.twiny()
            R, Nb = RudatToNbPhoton[0], RudatToNbPhoton[1]
            listX = []
            listRudat= [14,16,18,20,22,24]
            for rudat in listRudat:
                coeff = 10**((R-rudat)/10)
                listX.append(Nb*coeff)
            
            ax22.set_xlabel(r'Photon number, $n_{ph}$')
            newTickLocations = np.array(listRudat)
            ax22.set_xlim(ax2.get_xlim())
            ax22.set_xticks(newTickLocations)
            ax22.set_xticklabels(["{:.1f}".format(nb) for nb in listX])
        
        fig2.tight_layout()
        
        
        fig3 = plt.figure(figsize=(10,7))
        ax3 = fig3.add_subplot(1, 1, 1)
        ax3.grid() 
        
        h = 6.62607004*1e-34
        k = 1.38064852*1e-23
        f01 = 6.29911*1e9
        
        Ge_censor = np.array(dictJumpRate['excitation_rate_censor'])
        Gg_censor = np.array(dictJumpRate['decay_rate_censor'])
        
        Ge_proba = np.array(dictJumpRate['excitation_rate_proba'])
        Gg_proba = np.array(dictJumpRate['decay_rate_proba'])
        
        Teff_censor = - h*f01 / (k*np.log(Ge_censor/Gg_censor))
        Teff_proba = - h*f01 / (k*np.log(Ge_proba/Gg_proba))

        ax3.plot(self.listRudat, 1000*Teff_censor, marker = "P", ls = ':', color='#FF6C00', label='Exponential decaying')
        ax3.plot(self.listRudat, 1000*Teff_proba, marker = "s", ls = '-', color='green', label='Master equation')
       
        ax3.set_ylim(-400,450)
        ax3.set_xlabel('Rudat attenuation [dB]')
        ax3.set_ylabel('Effective temperature [mK]')
        ax3.legend(loc='upper right')
        
        if RudatToNbPhoton is not None:
            ax32 = ax3.twiny()
            R, Nb = RudatToNbPhoton[0], RudatToNbPhoton[1]
            listX = []
            listRudat= [14,16,18,20,22,24]
            for rudat in listRudat:
                coeff = 10**((R-rudat)/10)
                listX.append(Nb*coeff)
            
            ax32.set_xlabel(r'Photon number, $n_{ph}$')
            newTickLocations = np.array(listRudat)
            ax32.set_xlim(ax3.get_xlim())
            ax32.set_xticks(newTickLocations)
            ax32.set_xticklabels(["{:.1f}".format(nb) for nb in listX])
        
        fig3.tight_layout()

        if savepath != '':
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            plt.savefig(savepath+'JumpRates_vs_rudat_fit.png')
            plt.cla()
            plt.clf()
            plt.close(fig)
            gc.collect()
        else:
            plt.show()
            
