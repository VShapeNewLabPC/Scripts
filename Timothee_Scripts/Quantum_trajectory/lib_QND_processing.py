import os
import gc
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.optimize import minimize


from lib_Trajectory import Traject



T1 = 2000.0 #ns
MAX_BLANK = 40
MAX_BLANK_PLOT = MAX_BLANK-1
STATE_G = 0
STATE_E = 1

MEDIUM_SIZE = 18
BIGGER_SIZE = 20

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



def singleGaussian(x, amp, mean, std):
    return amp * np.exp(-1.0 * (x - mean)**2 / (2 * std**2))

def singleGaussianLabel(amp, mean, std):
    if mean < 0: 
        label = r'$%.3f exp(-(x + %.3f)^2 / (2*%.3f^2)$'%(amp,  abs(mean), std)
    else:
        label = r'$%.3f exp(-(x - %.3f)^2 / (2*%.3f^2)$'%(amp,  mean, std)
    return label
    

def doubleGaussian(x, amp1, mean1, std1, amp2, mean2, std2):
    return (singleGaussian(x, amp1, mean1, std1) + 
            singleGaussian(x, amp2, mean2, std2))

##############################################################
    
def exponential_decay(t, a, b):
    return a*np.exp(-b*t)

def exponential_decay_label(a, b):
    return r'$%.3fexp[-t/%.3f]$'%(a, 1.0/b)

##############################################################

def fitMeanTraject_func(t, t0, n0, G1, G2):
    G = G1 + G2
    egt = np.exp(-G*(t-t0))
    
    return n0*egt + G1*(1 - egt) / G 

def fitMeanTraject_funcLabel(t0, n0, G1, G2):
    return r'$n_o = %.3f$, $\Gamma'
    G = G1 + G2
    return r'$%.3fexp[-%.3e(t-%d)] + %.3f(1 - exp[-%.3e(t-%d)])/%.3f$'%(n0, G, t0, G1*1000, G, t0, G*1000)
    
# ##############################################################


def computeHistogram(data, nbins=None, normalized=True, naver=None, nbins_min=None, nbins_max=None):
    if nbins is None:
        fun_bins = 'auto'
        
        n = len(data)
        iqr = np.subtract(*np.percentile(data, [75, 25]))
        h = 2*iqr/n**(1/3)
        
        if h == 0:
            fun_bins = 5
        else:
            rangeData = np.max(data) - np.min(data)
            hbins = int(np.round(np.ceil(rangeData / h)))
            
            if naver is not None:
                if h < naver:
                    fun_bins = int(np.round(np.ceil(rangeData / naver)))
                    
            if nbins_max is not None:
                if hbins > nbins_max:
                    fun_bins = nbins_max
                    
            if nbins_min is not None:
                if hbins < nbins_min:
                    fun_bins = nbins_min
            
    else:
        fun_bins = nbins
   
    data_hist, bin_edges = np.histogram(data, bins=fun_bins)
    binscenters_hist = np.array([0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(len(bin_edges)-1)])
    if normalized:
        data_hist = 1.0*data_hist / len(data)

    return data_hist, binscenters_hist, bin_edges
    

def cumulativeFrequency(T, a, b, startT):
    out = []
    
    maxT = np.max(T)
    norma, error = quad(exponential_decay, startT, maxT, args=(a, b))
            
    t0 = T[0]
    result0, error = quad(exponential_decay, startT, t0, args=(a, b))
    
    for t in T:
        if t == t0:
            result = result0
        elif t == maxT:
            result = norma
            result0, t0 = result, t
        else:
            result, error = quad(exponential_decay, startT, t, args=(a, b)) 
            result0, t0 = result, t
           
        if norma == 0:
            out.append(0)
        else:
            out.append(result/norma)
        
    return out



############################################################################
###################   CLASS Base   #########################################
############################################################################

class Base(object):
    nameQNDfolder = None
    dataFolder = None
    rudat = None
    durationTraject = None
    
    listTypeTraject = None
    with_pio2 = None
    
    dictFileTraject = {'no': None, 'pi': None, 'pio2': None}
    dictAllTraject = {'no': None, 'pi': None, 'pio2': None}
    
    
    def __init__(self, dataFolder, with_pio2, tread, rudat):
        '''
        The initialisation function of the object Base

        Parameters
        ----------
        dataFolder : string
            the path to the folder containing the raw data folder.
        with_pio2 : bool
            if True, the measurements contain trajectories starting with a pi/2-pulse.
        tread : float
            the duration in microseconds [us] of a trajectory.
        rudat : float
            the value of rudat attenuation.

        Returns
        -------
        None.

        '''
        self._init_dataFolder(dataFolder)
            
        self.rudat = rudat
        
        self.with_pio2 = with_pio2
        self.durationTraject = tread*1000 #ns

        #initialisation
        self.listTypeTraject = ['no', 'pi']
        if self.with_pio2:
            self.listTypeTraject += ['pio2']

        self._init_dict()
     
        
    def _init_dict(self):
        for type_trj in self.listTypeTraject:
            self.dictFileTraject[type_trj] = []
            self.dictAllTraject[type_trj] = []
        
    def _init_dataFolder(self, dataFolder):
        if type(dataFolder) != str:
            self.dataFolder = [folder if folder.endswith('\\') 
                               else folder+'\\' for folder in dataFolder]
        else:
            self.dataFolder = [dataFolder if dataFolder.endswith('\\') 
                               else dataFolder+'\\']
                
        
    ################################
    ###   ADDING OF TRAJECTORY #####
    ################################
    
    def addTraject(self, traj):
        self.dictAllTraject[traj.type_trj].append(traj)


    def addTraject_fromData(self, type_trj, args_data, args_loading):
        traj = Traject(type_trj, args_data, args_loading)
        self.dictAllTraject[type_trj].append(traj)
        
    
    ####################################
    ###   SETTING ALL TRAJECTORIES #####
    ####################################
    
    def setTrajectMethod(self, name_method, comment=True, **kwargs):
        '''
        Take the name of a method defined in the object Traject 
        and apply it to all trajectories contained in self.dictAllTraject

        Parameters
        ----------
        name_method : string
            The name of the method to use.
        comment : Bool, optional
            If True, this function print a comment in the console.
            The default is True.
        **kwargs : method properties
            The argument that the method needs to use, cqn be empty.

        Returns
        -------
        None.

        '''
        nblank = MAX_BLANK
        
        if comment:
            start = time.time()
            sys.stdout.write('Set of {:{blank}} ...'.format(name_method, blank=nblank))
        
        for type_trj in self.listTypeTraject:
            for trj in self.dictAllTraject[type_trj]:
                getattr(trj, name_method)(**kwargs)
                
        if comment:
            end = time.time()
            sys.stdout.write('\rSet of {:{blank}} ...   done - {:.3f}s\n'.format(name_method, end-start, blank=nblank))


    def setTrajectAttribut(self, name_attr, val_attr):
        '''
        Take the name of the attribut defined in the object Traject
        and modify it using the given value

        Parameters
        ----------
        name_attr : string
            The name of the attribut to modify.
        val_attr : float
            The new value to apply to the attribut.

        Returns
        -------
        None.

        '''
        for type_trj in self.listTypeTraject:
            for trj in self.dictAllTraject[type_trj]:
                trj.__setattr__(name_attr, val_attr)
                
       
    def setNameQNDFolder(self, **kwargs):
        definition = self.get_args("definition", **kwargs)
        if definition is None:
            print("The argument 'definition' was not defined, default: 'Me'")
            definition = 'Me'
            
        threshold = self.get_args("threshold", **kwargs)
        if threshold is None:
            print("The argument 'threshold' was not defined, default: 'double'")
            threshold = 'double'
            
        security = self.get_args("security", **kwargs)
        if security is None:
            print("The argument 'security' was not defined, default: 'False'")
            security = 'False'
        
        self.nameQNDfolder = 'QNDness_'
        if threshold == 'double':
            self.nameQNDfolder += 'double_{}'.format(definition)
        else:
            self.nameQNDfolder += 'single_{}'.format(threshold)
            
        if not security:
            self.nameQNDfolder += '_bypass'
         
    #########################
    ###   GET FUNCTION  #####
    #########################
    
    def getOneBigList_of(self, str_attr, listType, startData=0, endData=None, addBefore_TypeTraject=False):
        big_list = np.array([])
        for type_trj in listType:
            for trj in self.dictAllTraject[type_trj]:
                listData = getattr(trj, str_attr)[startData:endData]
                
                if addBefore_TypeTraject:
                    if trj.type_trj == 'no':
                        listData = np.concatenate(([-1], listData))
                    if trj.type_trj == 'pi':
                        listData = np.concatenate(([-3], listData))
                    if trj.type_trj == 'pio2':
                        listData = np.concatenate(([-1.5], listData))
                    
                big_list = np.concatenate((big_list, listData))                                          
                
        return big_list
    
        
    def getNaverList(self):
        
        def getNaver_fromFile(file):
            end = file.partition('naver')[2]
            str_Naver = end.partition('.')[0]
            
            return int(str_Naver) 
        
        path = self.dataFolder[0] + 'Processed_data\\'
        try:
            listdir = os.listdir(path)
        except FileNotFoundError:
            print("There is no processed data, list_naver needs to be set manually")
            return []
            
        listAllNaver = sorted([getNaver_fromFile(file) for file in listdir])

        return listAllNaver
    
    def get_args(self, nameArgs, **kwargs):
        try:
            args = kwargs[nameArgs]
        except KeyError:
            return None
        
        return args
  
      
    
############################################################################
###################   CLASS Loading   ######################################
############################################################################

class Loading(Base):
    n_aver = None
    
    indStartTime = 0
    indEndTime = None
    
    loading_fromRawData = False
    loading_fromProcessedData = False
    
    
    ##################################
    ###   LOADING PROCEDURE ##########
    ##################################
    
    def loadingProcedure(self, n_aver):
        self.n_aver = n_aver
        if self.loading_fromRawData:
            self.setTrajectMethod('setNormalizedData', n_aver=n_aver, 
                                  ind_start=self.indStartTime, ind_stop=self.indEndTime)
        else:
            self._init_dict()
            loading = self.loadProcessedData()
            if not loading:
                self.loadRawData()
                self.setTrajectMethod('setNormalizedData', n_aver=n_aver, 
                                      ind_start=self.indStartTime, ind_stop=self.indEndTime)
            


    ####################################
    ###   LOADING OF RAW DATA ##########
    ####################################
    
    def loadRawData(self):
        print('Start to load raw trajectories')
                    
        self.loading_fromProcessedData = False
        self.loading_fromRawData = True
        self.getAdressRawData()
        for dtype in ['pi', 'no', 'pio2']:
            self.loadTraject_fromAdress(dtype)
        print(' ') 
        
        
    def getAdressRawData(self):

        def recognize_dir(dirname):
            '''
            function to check is a given name of directory correspond to parameters or data
            '''
            if 'traj_pio2' in dirname:
                return 'pio2'
            elif 'traj_no' in dirname:
                return 'no'
            elif 'traj_pi' in dirname:
                return 'pi'
            else:
                return 'nope'

        import os
        QND_data_folder = self.dataFolder[0] + 'Data\\'
            
        listdir = os.listdir(QND_data_folder)
        for fname in listdir:
            type_trj = recognize_dir(fname)
            if type_trj == 'nope':
                continue

            if fname.endswith('.dat'):
                self.dictFileTraject[type_trj].append( QND_data_folder + fname )

    
    def loadTraject_fromAdress(self, type_trj):

        if self.dictFileTraject[type_trj] == None:
            return False

        if len(self.dictFileTraject[type_trj]) == 0:
            return False

        n_max = len(self.dictFileTraject[type_trj])
        time_init = time.time()
        time_during = time_init
        sys.stdout.write('\r{} : {:03d}/{} ---- loading time : {:.2f}s ---- mean loading time : {:.2f}s'.format(type_trj, 
                             0, n_max, 0, 0))
        
        for n, file in enumerate(self.dictFileTraject[type_trj]):
            
            data = np.loadtxt(file)
            
            (STRINGS, ROWS) = np.shape(data)
            n_times = len(np.unique(data[:,0]))
            N_BLOCKS = int(STRINGS/n_times)
            data = np.reshape(data, (N_BLOCKS, n_times, ROWS))
            
            for num_trj in range(N_BLOCKS):
                args_data = [data[num_trj,:,0], data[num_trj,:,1], data[num_trj,:,2]]
                args_loading = [self.n_aver]
                self.addTraject_fromData(type_trj, args_data, args_loading)

                gc.collect()
            time_during = time.time()
            sys.stdout.write('\r{} : {:03d}/{} ---- loading time : {:.2f}s ---- mean loading time : {:.2f}s'.format(type_trj, 
                             n+1, n_max, time_during-time_init, (time_during-time_init)/(n+1)))
        print(' ')
        gc.collect()

    
    ##########################################
    ###   LOADING OF PROCESSED DATA ##########
    ##########################################
    
    def loadProcessedData(self):
        print('Start the process at {}'.format(str(datetime.now())))
        nblank = MAX_BLANK
        
        for folder in self.dataFolder:
            print('Folder of QND : {}'.format(folder))
            sys.stdout.write('\nTry to {:{blank}} ...'.format('load processed trajectories', blank=nblank))
            start = time.time()
            processed_data_path = self.searchProcessedData(folder)
            if processed_data_path:
                self.loading_fromRawData = False
                self.loading_fromProcessedData = True
                self.loadTraject_fromProcessedData(processed_data_path)
                end = time.time()
                sys.stdout.write('\rTry to {:{blank}} ...   done - {:.3f}s\n\n'.format('load processed trajectories', end-start, blank=nblank))
            else:
                sys.stdout.write('\rTry to load processed trajectories ... failed\n')
                return False
            
        return True
    
    def searchProcessedData(self, folder):
        import os
        processed_data_path = folder + 'Processed_data\\'
        if not os.path.exists(processed_data_path):
            return ''
          
        listdir = os.listdir(processed_data_path)
        for fname in listdir:
            if fname.endswith('naver{}.dat'.format(self.n_aver)):
                return processed_data_path + fname
            
        return ''
    
    
    def loadTraject_fromProcessedData(self, processed_data_path):        
        data = np.loadtxt(processed_data_path)
            
        (STRINGS, ROWS) = np.shape(data)
        unique_data = np.unique(data[:,0])
        if len(np.where(unique_data == -1.5)[0]) == 0:
            n_remove = 1
        else:
            n_remove = 2
        
        n_times = len(unique_data) - n_remove
        N_BLOCKS = int(STRINGS/n_times)
        data = np.reshape(data, (N_BLOCKS, n_times, ROWS))
        
        args_data = [None, None, None]
        args_loading = [self.n_aver]
        
        for num_trj in range(N_BLOCKS):
            time_comp = data[num_trj, :, 0]
            x_comp = data[num_trj, :, 1]
            if time_comp[0] == -1:
                type_trj = 'no'
            elif time_comp[0] == -3:
                type_trj = 'pi'
            else:
                type_trj = 'pio2'
            
            traj = Traject(type_trj, args_data, args_loading)
            traj.time_comp = time_comp[1:None]
            traj.x_comp = x_comp[1:None]
            self.addTraject(traj)
            
            
    def saveProcessedData(self, savepath, savename='Processed_data.dat'):
        if len(self.dataFolder) != 1:
            return
        
        if self.n_aver < 10:
            return
        
        #if data comes from a processed file, no need to recreate it
        if self.loading_fromProcessedData:
            return
        
        import os
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        
        listfile = os.listdir(savepath)
        if savename in listfile:
            return
        
        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Save of {:{blank}} ...'.format(savename, blank=nblank))
        
        list_time = self.getOneBigList_of('time_comp', self.listTypeTraject, addBefore_TypeTraject=True)
        list_x = self.getOneBigList_of('x_comp', self.listTypeTraject, addBefore_TypeTraject=True)
        
        mat_data = np.column_stack((list_time, list_x))
        
        header = 'time_comp [ns] \t x_comp [V]'
        np.savetxt(savepath+savename, mat_data, header=header)
        
        end = time.time()
        sys.stdout.write('\rSave of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
    
            
############################################################################
###################   CLASS MeanTraject   ##################################
############################################################################

class MeanTraject(Loading):
    dictStatTraject = None
    
    centerMean_g = None
    centerMean_e = None
    
    ind_centerMean_e = None
    ind_centerMean_g = None
    
    threshold_mean = None
    sign_ge = None
    
    coeffMean_pi = None
    coeffMean_no = None
    
    dictDurationStat = {'decay_rate_stat': None, 'excitation_rate_stat': None}

    
    def setMeanTraject(self):
        nblank = MAX_BLANK
        
        start = time.time()
        sys.stdout.write('Set of {:{blank}} ..'.format('statistics', blank=nblank))
        
        self.dictStatTraject = {'no': None, 'pi': None, 'pio2': None}
        for type_trj in self.listTypeTraject:
            self.dictStatTraject[type_trj] = {'x_mean': None, 'x_std': None}
            
            length_pt = len(self.dictAllTraject[type_trj][0].x_comp)
            length_trj = len(self.dictAllTraject[type_trj])
            
            if length_trj == 1:
                self.dictStatTraject[type_trj]['x_mean'] = self.dictAllTraject[type_trj][0].x_comp
                self.dictStatTraject[type_trj]['x_std'] = 0
                continue
                
            self.dictStatTraject[type_trj]['x_mean'] = np.zeros(length_pt)
            self.dictStatTraject[type_trj]['x_std'] = np.zeros(length_pt)

            for num_pt in range(length_pt):
                list_pt_trj = []
                for num_trj in range(length_trj):
                    list_pt_trj.append(self.dictAllTraject[type_trj][num_trj].x_comp[num_pt])

                self.dictStatTraject[type_trj]['x_mean'][num_pt] = np.mean(list_pt_trj)
                self.dictStatTraject[type_trj]['x_std'][num_pt] = np.std(list_pt_trj)
                
        end = time.time()
        sys.stdout.write('\rSet of {:{blank}} ...   done - {:.3f}s\n'.format('statistics', end-start, blank=nblank))

                    
    def setSingleThreshold(self):
        nblank = MAX_BLANK
        
        start = time.time()

        time_vec = self.dictAllTraject['pi'][0].time_comp
        time_min = 500
        condition = [False]
        while not np.any(condition):
            condition = np.where(time_vec<time_min, True, False)
            time_min = 2*time_min
           
        extractedMeanVector_e = self.dictStatTraject['pi']['x_mean'][condition]
        extractedMeanVector_g = self.dictStatTraject['no']['x_mean'][condition]
            
        mean_e = np.mean(extractedMeanVector_e)
        mean_g = np.mean(extractedMeanVector_g)
        
        if mean_e > mean_g:
            self.sign_ge = +1
            self.centerMean_e = np.max(extractedMeanVector_e)
            self.ind_centerMean_e = np.argmax(extractedMeanVector_e)
            self.centerMean_g = np.min(extractedMeanVector_g)
            self.ind_centerMean_g = np.argmin(extractedMeanVector_g)
        else:
            self.sign_ge = -1
            self.centerMean_e = np.min(extractedMeanVector_e)
            self.ind_centerMean_e = np.argmin(extractedMeanVector_e)
            self.centerMean_g = np.max(extractedMeanVector_g)
            self.ind_centerMean_g = np.argmax(extractedMeanVector_g)
            
        #Trajectory starts around 150ns and finish around tread+100ns
        if self.n_aver == 0:
            self.indStartTime = np.max((self.ind_centerMean_g, self.ind_centerMean_e))
            self.indEndTime = np.argwhere(time_vec == self.durationTraject)[0,0] + 50
            #50 due to the fact that trajectories does not start at t=0 but more around t=80ns
            
        self.threshold_mean = np.mean((self.centerMean_e, self.centerMean_g))
        
        end = time.time()
        sys.stdout.write('\r\rSet of {:{blank}} ...   done - {:.3f}s\n'.format('best threshold', end-start, blank=nblank))
        return True
        
    
    def fitMeanTraject(self):
        nblank = MAX_BLANK
        
        start = time.time()
        sys.stdout.write('Fit of {:{blank}} ...'.format('statistical {} trajectory'.format('all'), blank=nblank))

        if self.dictStatTraject['pi']['x_mean'] is None or self.dictStatTraject['no']['x_mean'] is None:
            self.setMeanTraject()
            
        time_vec = self.dictAllTraject['pi'][0].time_comp
        
        
        shift = np.mean([self.centerMean_e, self.centerMean_g])
        coeff = 1/abs(self.centerMean_e-self.centerMean_g)
        
        comboX = np.append(time_vec, time_vec)
        comboY = np.append((self.dictStatTraject['pi']['x_mean']-shift)*coeff + 0.5, 
                           (self.dictStatTraject['no']['x_mean']-shift)*coeff + 0.5)
        
        def FuncPen(comboData, n0, Gp, Gm, lenY1, t0):
            extract_e = comboData[0:lenY1] # first data
            extract_g = comboData[lenY1:None] # second data
            
            result_e = fitMeanTraject_func(extract_e, t0, n0, Gp, Gm)
            result_g = fitMeanTraject_func(extract_g, t0, 1-n0, Gp, Gm)
                        
            return np.append(result_e, result_g)
        
        funcFit = lambda comboData, n0, Gp, Gm: FuncPen(comboData, n0, Gp, Gm, len(time_vec), time_vec[0])
        
        coeff, error = curve_fit(funcFit, comboX, comboY, p0=(1, 1/T1, 1/T1), bounds=(0, +np.inf))
        n0, Gp, Gm = coeff
        
        self.dictDurationStat['decay_rate_stat'] = Gm
        self.dictDurationStat['excitation_rate_stat'] = Gp
        
        self.coeffMean_pi = (n0, Gp, Gm)
        self.coeffMean_no = (1-n0, Gp, Gm)
        
        end = time.time()
        sys.stdout.write('\rFit of {:{blank}} ...   done - {:.3f}s\n'.format('statistical {} trajectory'.format('all'), end-start, blank=nblank))
        
        
    def plotMeanTraject(self, save=True, savepath='', savename='Statistics_traj.png'):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        
        if self.dictStatTraject is None:
            sys.stdout.write('\rPlot of {:{blank}} ...   failed, no data\n'.format(savename, blank=nblank))
            return False
        
        if self.n_aver == 0:
            time_vec = self.dictAllTraject['no'][0].time
        else:
            time_vec = self.dictAllTraject['no'][0].time_comp

        fig = plt.figure(figsize=[18,8])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
        for type_trj, color_trj in zip(self.listTypeTraject, ['b', 'r', 'g']):
            ax.plot(time_vec, self.dictStatTraject[type_trj]['x_mean'], '-', color=color_trj, linewidth=1.5, label=type_trj)
            ax.fill_between(time_vec, self.dictStatTraject[type_trj]['x_mean']-self.dictStatTraject[type_trj]['x_std'],
                                      self.dictStatTraject[type_trj]['x_mean']+self.dictStatTraject[type_trj]['x_std'],
                                      color=color_trj, alpha=0.4)
        
        if self.n_aver == 0:
            ax.axvline(x=time_vec[self.indStartTime], ls=':', color='k')
            ax.axvline(x=time_vec[self.indEndTime], ls=':', color='k')
            
        spaceTime = np.linspace(time_vec[0], time_vec[-1], 5000)
            
        coeff = 1/abs(self.centerMean_e-self.centerMean_g)
        shift = np.mean([self.centerMean_e, self.centerMean_g])

        dataFit_pi = (fitMeanTraject_func(spaceTime, time_vec[0], *self.coeffMean_pi)-0.5)/coeff + shift
        dataFit_no = (fitMeanTraject_func(spaceTime, time_vec[0], *self.coeffMean_no)-0.5)/coeff + shift
        
        label_pi = r'$\Gamma^- = %.3f \mu s^{-1}$'%(self.coeffMean_pi[2]*1000)
        label_no = r'$\Gamma^+ = %.3f \mu s^{-1}$'%(self.coeffMean_pi[1]*1000)
        
        ax.plot(spaceTime, dataFit_no, '--b', label=label_no)
        ax.plot(spaceTime, dataFit_pi, '--r', label=label_pi)
            
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Mean normalized trajectorie [mV]')

        ax.legend()
        if savepath != '' and save:
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
            
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))


    def saveMeanTraject(self, savepath, savename='Statistics_traj.dat'):
        import os
        if not os.path.exists(savepath):
            os.makedirs(savepath)

        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Save of {:{blank}} ...'.format(savename, blank=nblank))

        header = ''
        length_pt = len(self.dictStatTraject['no']['x_mean'])
        mat_value = np.zeros((length_pt, 6))

        for type_trj in ['no', 'pi', 'pio2']:
            header += 'x_{}_mean'.format(type_trj) + '\t\t'
            header += 'x_{}_std'.format(type_trj) + '\t\t'

        for i, type_trj in enumerate(self.listTypeTraject):
            mat_value[:, 2*i] = self.dictStatTraject[type_trj]['x_mean']
            mat_value[:, 2*i+1] = self.dictStatTraject[type_trj]['x_std']

        np.savetxt(savepath+savename, mat_value, header=header)
        
        end = time.time()
        sys.stdout.write('\rSave of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))


############################################################################
###################   CLASS BigHistogram   #################################
############################################################################

class BigHistogram(MeanTraject):
    centerFit_g = None
    centerFit_e = None
    
    threshold_fit = None
    threshold_low = None
    threshold_high = None
    
    dataBigHistogram = None
    binsBigHistogram = None
    binsEdgesBigHistogram = None
    
    coeffGaussian_g = None
    coeffGaussian_e = None
    
    def setDoubleThreshold(self, **kwargs):
        definition = self.get_args("definition", **kwargs)
        if definition is None:
            print("The argument 'definition' was not defined, default: 'Me'")
            definition = 'Me'
        
        nblank = MAX_BLANK
        
        start = time.time()
        sys.stdout.write('Set of {:{blank}} ...'.format('improve threshold', blank=nblank))
        
        coeffLowGaussian, coeffHighGaussian = self.fitDoubleGaussian()
        
        amp_h, mean_h, std_h = coeffHighGaussian
        amp_l, mean_l, std_l = coeffLowGaussian
        
        # SNR = (mean_h - mean_l) / (std_h + std_l)
        #Slichter definition of high and low threshold
        if definition == 'Sl':
            self.threshold_fit = (mean_h + mean_l)/2
            self.threshold_high = self.threshold_fit + std_h**2 / (mean_h - mean_l)
            self.threshold_low = self.threshold_fit - std_l**2 / (mean_h - mean_l)
        
        #My definition of high and low threshold
        if definition == 'Me':
            self.threshold_fit = (mean_h + mean_l)/2
            self.threshold_high = (self.threshold_fit + mean_h) / 2
            self.threshold_low = (self.threshold_fit + mean_l) / 2
        
        if self.sign_ge > 0:
            self.coeffGaussian_g = coeffLowGaussian
            self.coeffGaussian_e = coeffHighGaussian
            self.centerFit_g, self.centerFit_e = mean_l, mean_h
        else:
            self.coeffGaussian_g = coeffHighGaussian
            self.coeffGaussian_e = coeffLowGaussian
            self.centerFit_g, self.centerFit_e = mean_h, mean_l
        
        end = time.time()
        sys.stdout.write('\rSet of {:{blank}} ...   done - {:.3f}s\n'.format('improve threshold', end-start, blank=nblank))
    
    def fitSingleTypeTraject(self, dtype, lowThreshold, highThreshold):
        sizeTraject = len(self.dictAllTraject[dtype][0].x_comp)
        bigTraject = self.getOneBigList_of('x_comp', [dtype], startData=0, endData=sizeTraject//2)
        
        dataBigHistogram, binsBigHistogram, binsEdgesBigHistogram = computeHistogram(bigTraject)
        
        lowCondition = np.where(binsBigHistogram<=lowThreshold, True, False)
        lowData = bigTraject[np.where(bigTraject<=lowThreshold, True, False)]
        lowGaussian = dataBigHistogram[lowCondition]
        
        highCondition = np.where(binsBigHistogram>=highThreshold, True, False)
        highData = bigTraject[np.where(bigTraject>=highThreshold, True, False)]
        highGaussian = dataBigHistogram[highCondition]        
        
        lowMeanMax = (self.threshold_mean + lowThreshold)/2
        highMeanMin = (self.threshold_mean + highThreshold)/2

        plow = [np.max(lowGaussian), np.mean(lowData), np.std(lowData)]
        phigh = [np.max(highGaussian), np.mean(highData), np.std(highData)]
        
        lowBoundsMin = [0, -np.inf, 0]
        lowBoundsMax = [+np.inf, lowMeanMax, +np.inf]
        
        highBoundsMin = [0, highMeanMin, 0]
        highBoundsMax = [+np.inf, +np.inf, +np.inf]
        
        guess, boundsMin, boundsMax = plow+phigh, lowBoundsMin+highBoundsMin, lowBoundsMax+highBoundsMax
        coeffDoubleGaussian, error =  curve_fit(doubleGaussian, binsBigHistogram, dataBigHistogram, 
                                                p0=guess, bounds=(boundsMin, boundsMax))
        
        return coeffDoubleGaussian[0:3], coeffDoubleGaussian[3:6]
       
        
    def fitDoubleGaussian(self):
                
        bigTraject = self.getOneBigList_of('x_comp', self.listTypeTraject)
        self.dataBigHistogram, self.binsBigHistogram, self.binsEdgesBigHistogram = computeHistogram(bigTraject)    
        
        if self.sign_ge > 0:
            lowThreshold = self.centerMean_g
            highThreshold = self.centerMean_e
            _, coeffHighGaussian = self.fitSingleTypeTraject('pi', lowThreshold, highThreshold)
            coeffLowGaussian, _  = self.fitSingleTypeTraject('no', lowThreshold, highThreshold)
        else:
            lowThreshold = self.centerMean_e
            highThreshold = self.centerMean_g
            coeffLowGaussian, _  = self.fitSingleTypeTraject('pi', lowThreshold, highThreshold)
            _, coeffHighGaussian = self.fitSingleTypeTraject('no', lowThreshold, highThreshold)
 
        lowBoundsMin = [0      , 1.05*coeffLowGaussian[1], 0.95*coeffLowGaussian[2]]
        lowBoundsMax = [+np.inf, 0.95*coeffLowGaussian[1], 1.05*coeffLowGaussian[2]]
        
        highBoundsMin = [0      , 0.95*coeffHighGaussian[1], 0.95*coeffHighGaussian[2]]
        highBoundsMax = [+np.inf, 1.05*coeffHighGaussian[1], 1.05*coeffHighGaussian[2]]
        
        guess = list(coeffLowGaussian) + list(coeffHighGaussian)
        boundsMin, boundsMax = lowBoundsMin+highBoundsMin, lowBoundsMax+highBoundsMax
        
        coeffDoubleGaussian, error =  curve_fit(doubleGaussian, self.binsBigHistogram, 
                                        self.dataBigHistogram, p0=guess)#, bounds=(boundsMin, boundsMax))
        
        coeffLowGaussian = coeffDoubleGaussian[0:3]
        coeffHighGaussian = coeffDoubleGaussian[3:6]
        
        return coeffLowGaussian, coeffHighGaussian

    
    def plotHistogramBigTraject(self, log=False, save=True, savepath='', savename='Hist_big_traj.png'):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Plot of {:{blank}} ...'.format(savename, blank=nblank))
        
        if self.dataBigHistogram is None or self.binsBigHistogram is None:
            sys.stdout.write('\rPlot of {:{blank}} ...   failed, no data\n'.format(savename, blank=nblank))
            return False
        
        bigTraject = self.getOneBigList_of('x_comp', self.listTypeTraject)
        mini, maxi = np.min(bigTraject), np.max(bigTraject)
        histSpace = np.linspace(mini, maxi, 5000)
            
        fig = plt.figure(figsize=[10, 7])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
        
        ax.bar(self.binsBigHistogram, self.dataBigHistogram, color='grey', log=log)

        ax.axvline(x=self.threshold_mean, ls=':', color='peru')
        ax.axvline(x=self.threshold_fit, color='green')
        ax.axvline(x=self.threshold_high, color='k',
                    label='High threshold : {:.3f}'.format(self.threshold_high))
        ax.axvline(x=self.threshold_low, color='k',
                    label='Low threshold : {:.3f}'.format(self.threshold_low))
        
        ax.axvline(x=self.centerMean_e, ls=':', color='r')
        ax.axvline(x=self.centerMean_g, ls=':', color='b')
        ax.axvline(x=self.centerFit_e, color='r')
        ax.axvline(x=self.centerFit_g, color='b')
        
        coeffDoubleGaussian = list(self.coeffGaussian_g) + list(self.coeffGaussian_e)
        ax.plot(histSpace, doubleGaussian(histSpace, *coeffDoubleGaussian), color='k')
        ax.plot(histSpace, singleGaussian(histSpace, *self.coeffGaussian_g), 
                color='b', label=singleGaussianLabel(*self.coeffGaussian_g))
        ax.plot(histSpace, singleGaussian(histSpace, *self.coeffGaussian_e), 
                color='r', label=singleGaussianLabel(*self.coeffGaussian_e))
        
        ax.set_xlim(-45,30)
        ax.set_xlabel('Voltage [mV]')
        ax.set_ylabel('Probability')
        ax.legend()
        
        if savepath != '' and save:
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
        
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
   
    
############################################################################
###################   CLASS StateDuration   ################################
############################################################################

class StateDuration(BigHistogram):
    
    dividerCensored_e = 20
    dividerAll_e = 20
    dividerCensored_g = 20
    dividerAll_g = 20
    
    listCensoredDuration_e = None
    listAllDuration_e = None
    listCensoredDuration_g = None
    listAllDuration_g = None
    
    histAllDuration_e = None
    histCensoredDuration_e = None
    histAllDuration_g = None
    histCensoredDuration_g = None
    
    binsAllDuration_e = None
    binsCensoredDuration_e = None
    binsAllDuration_g = None
    binsCensoredDuration_g = None
    
    coeffFitAllDuration_e = None
    coeffFitCensoredDuration_e = None
    coeffFitAllDuration_g = None
    coeffFitCensoredDuration_g = None
    
    dictDuration = {'decay_rate_All': None, 'excitation_rate_All': None,
                    'decay_rate_Censored': None, 'excitation_rate_Censored': None,
                    'meanDuration_e_all': None, 'meanDuration_g_all': None, 
                    'meanDuration_e_censored': None, 'meanDuration_g_censored': None}
    
    def setDurationList(self):

        # these lists contain censored durations (last duration is discarded)
        self.listCensoredDuration_e, self.listCensoredDuration_g = [], []
        # these lists contain all durations (no censoring)
        self.listAllDuration_e, self.listAllDuration_g = [], []
        
        for type_trj in self.listTypeTraject:
            for traj in self.dictAllTraject[type_trj]:                
                
                # this part is for cutting the tail of the trajectory, censoring
                ind_init = 0
                i = len(traj.ge_seq) - 1
                state_init = traj.ge_seq[i]
                while i>=0 and traj.ge_seq[i] == state_init:
                    i -= 1
                ind_final = i 
                #######
                
                #if trajectory contains only one state
                if ind_init > ind_final:
                    duration = traj.time_comp[-1] - traj.time_comp[0]
                    if traj.ge_seq[ind_init] == STATE_E:
                        self.listAllDuration_e.append(duration)
                    else:
                        self.listAllDuration_g.append(duration)
                    continue

                duration = traj.time_comp[-1] - traj.time_comp[ind_final]
                if traj.ge_seq[ind_init] == STATE_E:
                    self.listAllDuration_e.append(duration)
                else:
                    self.listAllDuration_g.append(duration)
                
                time_init, state_init = traj.time_comp[ind_init], traj.ge_seq[ind_init]
                for i in range(ind_init, ind_final):
                    state = traj.ge_seq[i]
                    
                    #if a jump occurs
                    if state != state_init:
                        
                       #       ________
                       #      /        \
                       # ____/<-------->\_______
                       #      time jump
                        
                        time_jump = (traj.time_comp[i-1] + traj.time_comp[i])/2
                        duration = time_jump - time_init
                        if state_init == STATE_E:
                            self.listCensoredDuration_e.append(duration)
                            self.listAllDuration_e.append(duration)
                        if state_init == STATE_G:
                            self.listCensoredDuration_g.append(duration)
                            self.listAllDuration_g.append(duration)

                        time_init, state_init = time_jump, state
                        
        
        self.dictDuration['meanDuration_e_censored'] = np.mean(self.listCensoredDuration_e)
        self.dictDuration['meanDuration_g_censored'] = np.mean(self.listCensoredDuration_g)
        self.dictDuration['meanDuration_e_all'] = np.mean(self.listAllDuration_e)
        self.dictDuration['meanDuration_g_all'] = np.mean(self.listAllDuration_g)

        
    def computeJumpRate(self):
        nblank = MAX_BLANK
        startTime = time.time()
        sys.stdout.write('Set of {:{blank}} ...'.format('Jump rates', blank=nblank))
        
        for state, rate in zip(['e', 'g'], ['decay_rate', 'excitation_rate']):
            
            for dtype in ['Censored', 'All']:
                listDuration = getattr(self, 'list{}Duration_{}'.format(dtype, state))
                divider = getattr(self, 'divider{}_{}'.format(dtype, state))

                # Fit of cumulative frequency
                #############################
                
                X = np.sort(listDuration)
                
                #cut the first part because not interesting due to false jumps
                mini = X[-1]/divider
                while mini < X[0]:
                    divider = divider - 3
                    mini = X[-1]/divider
                    if divider < 0:
                        divider = getattr(self, 'divider{}_{}'.format(dtype, state))
                        break
                    
                mini = X[-1]/divider
                
                condition = np.where(X>mini, True, False)
                Xfit = X[condition]
                Yfit = np.arange(1, len(Xfit)+1) / len(Xfit)
                
                function = lambda x, a, b: cumulativeFrequency(x, a, b, Xfit[0])
                coeffFitDuration, error = curve_fit(function, Xfit, Yfit, 
                                                    p0=(1, 1/T1), bounds=(0, +np.inf))
                    
                # Creation of histogram
                #######################
                
                if dtype == 'Censored':
                    histDuration, binsDuration, binsEdge = computeHistogram(listDuration, naver=self.n_aver, nbins_min=3, nbins_max=50)
                if dtype == 'All':
                    histDuration, binsDuration, binsEdge = computeHistogram(listDuration, nbins=binsEdge)
                
                # Fit of histogram
                ###################
                
                condition = np.where(histDuration!=0, True, False)
                histFit = histDuration[condition]
                binsFit = binsDuration[condition]
                
                if len(histFit) > 2:
                    start, stop = 1, -1
                else:
                    start, stop = 0, None
                    
                A0 = np.max(histFit[start:stop])   
                try:
                    function = lambda x, a : exponential_decay(x, a, coeffFitDuration[1])
                    a, error = curve_fit(function, binsFit[start:stop], histFit[start:stop], 
                                         p0=(A0), bounds=(0, +np.inf))
                except RuntimeError:
                    print('Curve fit of {} {} state duration histogram failed ...'.format(state, dtype))
                    a = 0
                
                coeffFitDuration[0] = a
                
                # Saving of result
                ##################
                
                self.dictDuration[rate+'_'+dtype] = coeffFitDuration[1] 
                setattr(self, 'hist{}Duration_{}'.format(dtype, state), histDuration)
                setattr(self, 'bins{}Duration_{}'.format(dtype, state), binsDuration)
                setattr(self, 'coeffFit{}Duration_{}'.format(dtype, state), coeffFitDuration)
                setattr(self, 'divider{}_{}'.format(dtype, state), divider)
                
        endTime = time.time()
        sys.stdout.write('\rSet of {:{blank}} ...   done - {:.3f}s\n'.format('Jump rates', endTime-startTime, blank=nblank))
    
        
    def plotHistogramDuration(self, fit=True, log=False, save=True, savepath='', savename='Hist_dtj.png'):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Plot of {:{blank}} ...'.format(savename, blank=nblank))
        
        if log:
            fig = plt.figure('Hist_dtj_log.png', figsize=[10,8])
        else:
            fig = plt.figure('Hist_dtj.png', figsize=[10,6])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid(b=True)
        
        dictColor = {'e': 'orange','g': 'cyan'}
        dictColorFit = {'e': 'red','g': 'blue'}
        
        dtype = 'Censored'
        for state in ['g', 'e']:
            histDuration = getattr(self, 'hist{}Duration_{}'.format(dtype, state))
            binsDuration = getattr(self, 'bins{}Duration_{}'.format(dtype, state))
            if histDuration is None or binsDuration is None:
                sys.stdout.write('\rPlot of {:{blank}} ...   failed, no data\n'.format(savename, blank=nblank))
                return False
            
            color = dictColor[state]
            ax.step([0]+list(binsDuration)+[0], [0]+list(histDuration)+[0], color=color, label=r'$|%s\rangle$ state'%(state))
            if fit:
                colorFit = dictColorFit[state]
                
                mini, maxi = np.min(binsDuration), np.max(binsDuration)
                hist_space = np.linspace(mini, maxi, 5000)
                
                coeffExp = getattr(self, 'coeffFit{}Duration_{}'.format(dtype, state))
                ax.plot(hist_space, exponential_decay(hist_space, *coeffExp),
                        color=colorFit)
                
        if log:
            ax.set_yscale("log")
        ax.set_xlabel('$\Delta t_{jump}$ [ns]')
        ax.set_ylabel('Probability')
        ax.legend()
        
        fig.tight_layout()
        
        if savepath != '' and save:
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
        
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
    
    
    def plotCumulativeFrequency(self, save=True, savepath='', savename='Cumulative_frequency.png'):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Plot of {:{blank}} ...'.format(savename, blank=nblank))
        
        dictColor = {'e': 'orange','g': 'cyan'}
        dictColorFit = {'e': 'red','g': 'blue'}
          
        fig = plt.figure('Cumulative_frequency', figsize=[10,6])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid(b=True)
        
        for state in ['e', 'g']:
        
            listCensoredDuration = getattr(self, 'listCensoredDuration_{}'.format(state))
            if listCensoredDuration is None:
                sys.stdout.write('\rPlot of {:{blank}} ...   failed, no data\n'.format(savename, blank=nblank))
                return False
            
            dividerCensored = getattr(self, 'dividerCensored_{}'.format(state))
            
            Xcensored = np.sort(listCensoredDuration)
            
            color = dictColor[state]
            colorFit = dictColorFit[state]
                        
            miniCensored = Xcensored[-1]/dividerCensored
            condition = np.where(Xcensored>miniCensored, True, False)
            
            XCensoredfit = Xcensored[condition]
            YCensoredfit = np.arange(1, len(XCensoredfit)+1) / len(XCensoredfit)
            
            ax.plot(XCensoredfit, YCensoredfit, color=color, label=r'$|%s\rangle$ state'%(state))
            
            coeffExpCensored = getattr(self, 'coeffFitCensoredDuration_{}'.format(state))
            coeff = list(coeffExpCensored) + [miniCensored]
            XCensoredSpace = np.linspace(np.min(XCensoredfit), np.max(XCensoredfit), 5000)
            ax.plot(XCensoredSpace, cumulativeFrequency(XCensoredSpace, *coeff),
                     '--', color=colorFit)

            
        ax.set_xlabel('$\Delta t_{jump}$ [ns]')
        ax.set_ylabel('Cumulative frequency')
        ax.legend()
        
        fig.tight_layout()
        
        if savepath != '' and save:
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
        
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
    
      
    def saveHistogramDuration(self, list_hist, savepath, savename='Hist_dtj.dat'):
        import os
        if not os.path.exists(savepath):
            os.makedirs(savepath)
            
        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Save of {:{blank}} ...'.format(savename, blank=nblank))
        
        header = 'duration [ns]'
        np.savetxt(savepath+savename, list_hist, header=header)
        
        end = time.time()
        sys.stdout.write('\rSave of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
  
    
############################################################################
###################   CLASS StateVoltage   #################################
############################################################################
          
class StateVoltage(StateDuration):
    stdVoltage_e = None
    stdVoltage_g = None
    
    histVoltage_e = None
    histVoltage_g = None
    
    binsVoltage_e = None
    binsVoltage_g = None
    
    coeffVoltage_e = None
    coeffVoltage_g = None
    
    errorVoltage_e = None
    errorVoltage_g = None
    
    def setHistogram_ge(self):
        nblank = MAX_BLANK
        
        start = time.time()
        sys.stdout.write('Set of {:{blank}} ...'.format('e-g histogram', blank=nblank))
        
        bigTraject = self.getOneBigList_of('x_comp', self.listTypeTraject)
        big_ge_seq = self.getOneBigList_of('ge_seq', self.listTypeTraject)
        
        condition_e = np.where(big_ge_seq==STATE_E, True, False)
        condition_g = np.logical_not(condition_e)
        
        listVoltage_e = np.extract(condition_e, bigTraject)
        listVoltage_g = np.extract(condition_g, bigTraject)
        
        self.stdVoltage_e = np.std(listVoltage_e)
        self.stdVoltage_g = np.std(listVoltage_g)
                
        self.histVoltage_e, self.binsVoltage_e, _ = computeHistogram(listVoltage_e)
        self.histVoltage_g, self.binsVoltage_g, _ = computeHistogram(listVoltage_g)
        
        p0 = [np.max(self.histVoltage_e), np.mean(listVoltage_e), np.std(listVoltage_e)]
        self.coeffVoltage_e, error =  curve_fit(singleGaussian, self.binsVoltage_e, self.histVoltage_e, p0=p0)
        self.errorVoltage_e = np.sqrt(np.diag(error))
        
        p0 = [np.max(self.histVoltage_g), np.mean(listVoltage_g), np.std(listVoltage_g)]
        self.coeffVoltage_g, error =  curve_fit(singleGaussian, self.binsVoltage_g, self.histVoltage_g, p0=p0)
        self.errorVoltage_g = np.sqrt(np.diag(error))
        
        end = time.time()
        sys.stdout.write('\rSet of {:{blank}} ...   done - {:.3f}s\n'.format('e-g histogram', end-start, blank=nblank))
    

    def plotHistogram_ge(self, log=False, save=True, savepath='', savename='Hist_eg_state.png'):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Plot of {:{blank}} ...'.format(savename, blank=nblank))
                    
        if self.histVoltage_g is None or self.histVoltage_e is None:
            sys.stdout.write('\rPlot of {:{blank}} ...   failed, no data\n'.format(savename, blank=nblank))
            return False
        
        fig = plt.figure(figsize=[18,8])

        for state, color, subplot in zip(['g', 'e'], ['b', 'r'], [121, 122]):
            binsVoltage = getattr(self, 'binsVoltage_{}'.format(state))
            histVoltage = getattr(self, 'histVoltage_{}'.format(state))
            coeffVoltage = getattr(self, 'coeffVoltage_{}'.format(state))
            errorVoltage = getattr(self, 'errorVoltage_{}'.format(state))
            
            ax = fig.add_subplot(subplot)    
            ax.grid()
        
            mini, maxi = np.min(binsVoltage), np.max(binsVoltage)
            histSpace = np.linspace(mini, maxi, 5000)
            ax.bar(binsVoltage, histVoltage, log=log, color=color)
            ax.plot(histSpace, singleGaussian(histSpace, *coeffVoltage), color='k')
            
            ax.set_xlabel('Voltage [mV]')
            ax.set_ylabel('Probability')
            ax.set_title('Histogram for $|{}\\rangle$ state'.format(state))
            
            label = r'$A=%.3f, A_{error}=%.2e$'%(coeffVoltage[0], errorVoltage[0]) +'\n'
            label += r'$\mu=%.3f, \mu_{error}=%.2e$'%(coeffVoltage[1], errorVoltage[1]) +'\n'
            label += r'$\sigma=%.3f, \sigma_{error}=%.2e$'%(coeffVoltage[2], errorVoltage[2])
            ax.text(0.7, 0.87, label, bbox=dict(facecolor='white'), transform=ax.transAxes)
            
        if savepath != '' and save:
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
        
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))      
   
  
############################################################################
###################   CLASS MeanProbaState   ###############################
############################################################################      
  
class MeanProbaState(StateVoltage):
    
    coeffProba_e = None
    coeffProba_g = None
    
    meanProba_e = None
    meanProba_g = None
    
    dictDurationProba = {'decay_rate_proba': None, 'excitation_rate_proba': None}
    
    def setMeanProba(self):
        nblank = MAX_BLANK
        start = time.time()
        sys.stdout.write('set of {:{blank}} ... '.format('mean probabilities', blank=nblank))
        
        length_pt = len(self.dictAllTraject['pi'][0].ge_seq)

        self.meanProba_e = np.zeros(length_pt)
        self.meanProba_g = np.zeros(length_pt)
        
        for num_pt in range(length_pt):
            nb_e = 0
            nb_g = 0
            for type_trj in self.listTypeTraject:
                nbTraject = len(self.dictAllTraject[type_trj])
 
                for num_trj in range(nbTraject):
                    if self.dictAllTraject[type_trj][num_trj].ge_seq[num_pt] == STATE_E:
                        nb_e += 1
                    else:
                        nb_g += 1
    
                self.meanProba_e[num_pt] = nb_e / (nb_e + nb_g)
                self.meanProba_g[num_pt] = nb_g / (nb_e + nb_g)
    
        end = time.time()
        sys.stdout.write('\rset of {:{blank}} ...   done - {:.3f}s\n'.format('mean probabilities', end-start, blank=nblank))
                
        
    
    def fitMeanProba(self):
        nblank = MAX_BLANK
        
        start = time.time()
        sys.stdout.write('Fit of {:{blank}} ...'.format('mean probabilities', blank=nblank))

        if self.meanProba_e is None or self.meanProba_g is None:
            self.setMeanProba()
            
        time_vec = self.dictAllTraject['pi'][0].time_comp
        
        comboX = np.append(time_vec, time_vec)
        comboY = np.append(self.meanProba_e, self.meanProba_g)
        
        def FuncPen(comboData, n0, Gp, Gm, lenY1, t0):
            extract_e = comboData[0:lenY1] # first data
            extract_g = comboData[lenY1:None] # second data
            
            result_e = fitMeanTraject_func(extract_e, t0, n0, Gp, Gm)
            result_g = fitMeanTraject_func(extract_g, t0, 1-n0, Gm, Gp)
            
            expectedResultSum = np.ones(np.shape(result_e))
            penalization = np.sum(expectedResultSum - (result_e+result_g))
            
            # return np.append(result_e, result_g)
            return np.append(result_e+penalization, result_g+penalization)
        
        funcFit = lambda comboData, n0, Gp, Gm: FuncPen(comboData, n0, Gp, Gm, len(time_vec), time_vec[0])
        
        coeff, error = curve_fit(funcFit, comboX, comboY, p0=(0.5, 1/T1, 1/T1), bounds=(0, +np.inf))
        n0, Gp, Gm = coeff
        
        self.dictDurationProba['decay_rate_proba'] = Gm
        self.dictDurationProba['excitation_rate_proba'] = Gp
        
        self.coeffProba_e = (n0, Gp, Gm)
        self.coeffProba_g = (1-n0, Gm, Gp)
        
        end = time.time()
        sys.stdout.write('\rFit of {:{blank}} ...   done - {:.3f}s\n'.format('mean probabilities', end-start, blank=nblank))
        
        
            
    def plotMeanProba(self, save=True, savepath='', savename='Mean_Probabilities.png'):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        
        fig = plt.figure(figsize=[10,6])
        ax = fig.add_subplot(1, 1, 1)
        ax.grid()
        
        time_vec = self.dictAllTraject['pi'][0].time_comp
        
        ax.plot(time_vec, self.meanProba_e, '-r', label=r'$|e\rangle$')
        ax.plot(time_vec, self.meanProba_g, '-b', label=r'$|g\rangle$')
        
        spaceTime = np.linspace(time_vec[0], time_vec[-1], 5000)
        fitData_e = fitMeanTraject_func(spaceTime, time_vec[0], *self.coeffProba_e)
        fitData_g = fitMeanTraject_func(spaceTime, time_vec[0], *self.coeffProba_g)
        
        label_e = r'$n_0 = %.3f$, $\Gamma^+ = %.3f \mu s^{-1}$'%(self.coeffProba_e[0], self.coeffProba_e[1]*1000)
        label_g = r'$n_0 = %.3f$, $\Gamma^- = %.3f \mu s^{-1}$'%(self.coeffProba_g[0], self.coeffProba_g[1]*1000)
        ax.plot(spaceTime, fitData_e, '--r', label=label_e)
        ax.plot(spaceTime, fitData_g, '--b', label=label_g)
        
        ax.legend()
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Probability')
        
        title = r'Probability to be in $|e\rangle$ and $|g\rangle$ states as function of time, rudat={}dB'.format(self.rudat)
        ax.set_title(title)
        
        fig.tight_layout()
        if savepath != '' and save:
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
            
        end = time.time()
        sys.stdout.write('Plot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
                
     
############################################################################
###################   CLASS Unknown   ######################################
############################################################################

class Distance_eg(MeanProbaState):
    
    dictMatData = None
    
    def computeDistance_eg(self):
        nblank = MAX_BLANK
        
        start = time.time()
        sys.stdout.write('Set of {:{blank}} ...'.format('unknown hist', blank=nblank))
        
        stopNaver = self.indEndTime - self.indStartTime
        listNaver = np.arange(1, stopNaver+1, 1)
        
        keys = self.listTypeTraject
        valueX = [[] for _ in self.listTypeTraject]
        valueY = [[] for _ in self.listTypeTraject]

        self.dictUnknown_Y = dict(zip(keys, valueY))
        self.dictUnknown_X = dict(zip(keys, valueX))
        for Naver in listNaver:
            self.setTrajectMethod('setNormalizedData', comment=False, n_aver=Naver, 
                                  ind_start=self.indStartTime, ind_stop=self.indStartTime+Naver)
            
            for type_trj in self.listTypeTraject:
                
                for trj in self.dictAllTraject[type_trj]:
                    
                    firstPoint = trj.x_comp[0]
                    self.dictUnknown_Y[type_trj].append(firstPoint)
                    self.dictUnknown_X[type_trj].append(Naver)
                
        
        self.setTrajectMethod('setNormalizedData', comment=False, n_aver=self.n_aver, 
                                  ind_start=self.indStartTime, ind_stop=self.indEndTime)
        
        end = time.time()
        sys.stdout.write('\rSet of {:{blank}} ...   done - {:.3f}s\n'.format('unknown hist', end-start, blank=nblank))
        
        
                
    def plotDistance_eg(self, save=True, savepath='', savename="Unknow plot"):
        nblank = MAX_BLANK_PLOT
        start = time.time()
        
        sys.stdout.write('Plot of {:{blank}} ...'.format(savename, blank=nblank))

        stopNaver = self.indEndTime - self.indStartTime
        listNaver = np.arange(1, stopNaver+1, 1)
        
        list_X = self.dictUnknown_X['no'] + self.dictUnknown_X['pi'] + self.dictUnknown_X['pio2']
        list_Y = self.dictUnknown_Y['no'] + self.dictUnknown_Y['pi'] + self.dictUnknown_Y['pio2']

        meanVoltage = np.mean(list_Y)
        stdVoltage = np.std(list_Y)
        
        cmap = "jet"
        
        fig, axes = plt.subplots(nrows=2, ncols=2, sharex='col',
                                 sharey='row', figsize=[18,8])
        
        h = axes[0,0].hist2d(self.dictUnknown_X['no'], self.dictUnknown_Y['no'], bins=len(listNaver), cmap=cmap)
        h = axes[0,1].hist2d(self.dictUnknown_X['pi'], self.dictUnknown_Y['pi'], bins=len(listNaver), cmap=cmap)
        h = axes[1,0].hist2d(self.dictUnknown_X['pio2'], self.dictUnknown_Y['pio2'], bins=len(listNaver), cmap=cmap)
        h = axes[1,1].hist2d(list_X, list_Y, bins=len(listNaver), cmap=cmap)

        for ax in axes.flat:
            ax.set(xlabel='Integration Time [ns]', ylabel='Voltage [mV]')
            ax.set_ylim(meanVoltage-3*stdVoltage, meanVoltage+3*stdVoltage)

        for ax in axes.flat:
            ax.label_outer()
        
        fig.colorbar(h[3], ax=axes.ravel().tolist(), orientation="vertical")
        
        if savepath != '' and save:
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
            
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
                
     
############################################################################
###################   CLASS QNDness   ######################################
############################################################################

class QNDness(Distance_eg):

    dictJump = {'n_jump': None, 'n_not_jump': None}
    dictQND = {'Q': None, 'Qe': None, 'Qg': None, 'uncertainty': None,
                'Pgg': None, 'Pee': None, 'Peg': None, 'Pge': None}
    
    def computeQNDness(self):
        """dictQND = {'Q': None, 'Qe': None, 'Qg': None,
                    'Pgg': None, 'Pee': None, 'Peg': None, 'Pge': None,
                    'n_jump': None, 'uncertainty': None}"""
    
        self.setTrajectMethod('setJumpCounter')
    
        n_jump, n_not_jump = 0, 0
        Pgg, Pee, Pge, Peg = 0, 0, 0, 0
        for trj in self.dictAllTraject['pi'] + self.dictAllTraject['no']:
            Pgg += trj.p_gg
            Pee += trj.p_ee
            Pge += trj.p_ge
            Peg += trj.p_eg
            
            n_jump += trj.n_jump
            n_not_jump += trj.n_not_jump

        
        self.dictJump['n_jump'] = n_jump
        self.dictJump['n_not_jump'] = n_not_jump

        self.dictQND['Pgg'] = float(Pgg)
        self.dictQND['Pee'] = float(Pee)
        self.dictQND['Pge'] = float(Pge)
        self.dictQND['Peg'] = float(Peg)

        Pg = self.dictQND['Pge'] + self.dictQND['Pgg']  #to normalize into probability
        Pe = self.dictQND['Pee'] + self.dictQND['Peg']  #to normalize into probability
        self.dictQND['Qg'] = self.dictQND['Pgg'] / Pg
        self.dictQND['Qe'] = self.dictQND['Pee'] / Pe
        self.dictQND['Q'] = (self.dictQND['Qg'] + self.dictQND['Qe'])/2.

        # we have an uncertainty on Q due to finite counting number of
        self.dictQND['uncertainty'] = (100./np.sqrt(Pg) + 100./np.sqrt(Pe))/2 # in percent
    
    
    def colorplotAllTraject(self, type_trj, vmin=None, vmax=None, save=True, savepath='', savename='Colorplot.png'):
        '''
        Takes a list of trajectories
        Plots two graphics:
        1) all trajectories. Voltage is color
        2) extracted states
        '''
        import matplotlib.gridspec as gridspec
        nblank = MAX_BLANK_PLOT
        
        start = time.time()
        sys.stdout.write('Plot of {:{blank}} ...'.format(savename, blank=nblank))
        
        list_trj = self.dictAllTraject[type_trj]
        time_list = list_trj[0].time_comp
        indt_list = range(len(list_trj))
        mV_2D_list = []
        state_2D_list = []

        for i in indt_list:
            trj = list_trj[i]
            mV_2D_list.append( trj.x_comp )
            state_2D_list.append( trj.ge_seq )


        #################################################
        fig = plt.figure(figsize=[20,40])
        gs = gridspec.GridSpec( 1, 2, width_ratios=[1,1])
        ax_raw = plt.subplot(gs[0])
        ax_state = plt.subplot(gs[1])

        c1 = ax_raw.pcolormesh(time_list, indt_list, mV_2D_list, vmin=vmin, vmax=vmax)
        c2 = ax_state.pcolormesh(time_list, indt_list, state_2D_list)
        fig.colorbar(c1, ax=ax_raw, orientation = 'horizontal')
        fig.colorbar(c2, ax=ax_state, orientation = 'horizontal')

        ax_raw.set_title('Normalized trajectories [mV]')
        ax_raw.set_ylabel('# of trajectorie')
        ax_raw.set_xlabel('Time [ns]')

        ax_state.set_title('State of qubit (g or e)')
        ax_state.set_xlabel('Time [ns]')

        ##################################################
        if savepath != '' and save:
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
            
        end = time.time()
        sys.stdout.write('\rPlot of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))

        ################
        return True
    
    
    def plotAllTraject(self, percent, savepath='', type_name=''):
        percent = percent/100 #to get the alue between 0 and 1
        
        if percent == 0:
            print('no plot to do')
            return True
        
        print(' ')
        for type_trj in self.listTypeTraject:
            
            start = time.time()
            n_len = len(self.dictAllTraject[type_trj])
            n_max = int(percent*n_len)
            for n in range(n_max):
                num_traj = int(n/percent)
                trj = self.dictAllTraject[type_trj][num_traj]
                sys.stdout.write('\r{} : plot saved {}/{}'.format(type_trj, n+1, n_max))
                traj_name = 'traj_{}_{}_{}'.format(type_trj, n+1, type_name)
                trj.plotter(savepath=savepath, savename=traj_name+'.png')
                
            end = time.time()
            sys.stdout.write('\rPlot of {} trajectories ...   done - {:.3f}s'.format(type_trj, end-start))
            print(' ')


    ###########################
    ###   SAVE Files  ######### 
    ###########################
    
    def saveQNDness(self, save, savepath, savename='QNDness.dat'):
        if not save:
            return False
        
        import os
        if not os.path.exists(savepath):
            os.makedirs(savepath)

        nblank = MAX_BLANK_PLOT
        start = time.time()
        sys.stdout.write('Save of {:{blank}} ...'.format(savename, blank=nblank))

        header = ''
        list_value = []
        listDict = [self.dictQND, self.dictJump, self.dictDuration, 
                    self.dictDurationProba, self.dictDurationStat]
        
        for ddict in listDict:
            for key, value in ddict.items():
                header += key + '\n'
                list_value.append(value)
            
        np.savetxt(savepath+savename, list_value, header=header)
        end = time.time()
        sys.stdout.write('\rSave of {:{blank}} ...   done - {:.3f}s\n'.format(savename, end-start, blank=nblank))
        

    #######################
    ###   DO ALL  #########
    #######################
    
    def setAllTraject_StateSeq(self, **kwargs):
        eg_center = self.get_args("eg_center", **kwargs)
        if eg_center is None:
            print("The argument 'eg_center' was not defined, default: 'fit'")
            eg_center = 'fit'
            
        threshold = self.get_args("threshold", **kwargs)
        if threshold is None:
            print("The argument 'threshold' was not defined, default: 'double'")
            threshold = 'double'
        
        size = self.get_args("size", **kwargs)
        if size is None:
            print("The argument 'size' was not defined, default: 'True'")
            size = True
            
        security = self.get_args("security", **kwargs)
        if security is None:
            print("The argument 'security' was not defined, default: 'False'")
            security = False
            
        
        if eg_center=='mean':
            self.setTrajectAttribut('e_center', self.centerMean_e)
            self.setTrajectAttribut('g_center', self.centerMean_g)
        else:
            self.setTrajectAttribut('e_center', self.centerFit_e)
            self.setTrajectAttribut('g_center', self.centerFit_g)
            
            
        if threshold == 'mean':
            self.setTrajectAttribut('threshold', self.threshold_mean)
            self.setTrajectMethod('setStateSeq_singleThreshold', 
                                  sign_ge=self.sign_ge, 
                                  threshold=self.threshold_mean,
                                  security=security)
        elif threshold == 'fit':
            self.setTrajectAttribut('threshold', self.threshold_fit)
            self.setTrajectMethod('setStateSeq_singleThreshold',
                                  sign_ge=self.sign_ge,
                                  threshold=self.threshold_fit,
                                  security=security)
        else:
            self.setTrajectAttribut('threshold_low', self.threshold_low)
            self.setTrajectAttribut('threshold_high', self.threshold_high)
            self.setTrajectMethod('setStateSeq_doubleThreshold', 
                                  sign_ge=self.sign_ge, 
                                  threshold_fit=self.threshold_fit, 
                                  threshold_low=self.threshold_low, 
                                  threshold_high=self.threshold_high,
                                  security=security)
         
            
        if size:
            self.setTrajectAttribut('e_size', self.stdVoltage_e)
            self.setTrajectAttribut('g_size', self.stdVoltage_g)
            
        
    def compute_all(self, **kwargs):
        saveOnlyProcessed = self.get_args("saveOnlyProcessed", **kwargs)
        if saveOnlyProcessed is None:
            print("The argument 'saveOnlyProcessed' was not defined, default: 'False'")
            saveOnlyProcessed = False
        
        self.setNameQNDFolder(**kwargs)
        
        if saveOnlyProcessed:
            self.setMeanTraject()
            self.setSingleThreshold()
            return True
        
        self.setMeanTraject()
        self.setSingleThreshold()
        self.fitMeanTraject()
        self.setDoubleThreshold(**kwargs) 
        self.setAllTraject_StateSeq(**kwargs)
        
        self.setDurationList()
        self.computeJumpRate()
        
        self.setHistogram_ge()
        
        self.setMeanProba()
        self.fitMeanProba()
        self.computeQNDness()
        
        if self.n_aver == 0:
            self.computeDistance_eg()
        
        
    def save_all(self, percent_plot=0.1, **kwargs):
        saveOnlyProcessed = self.get_args("saveOnlyProcessed", **kwargs)
        if saveOnlyProcessed is None:
            print("The argument 'saveOnlyProcessed' was not defined, default: 'False'")
            saveOnlyProcessed = False
        
        saveResult = self.get_args("saveResult", **kwargs)
        if saveResult is None:
            print("The argument 'saveResult' was not defined, default: 'True'")
            saveResult = True
        
        if len(self.dataFolder) == 1:
            savepath = self.dataFolder[0]
        else:
            listFolder = self.dataFolder[0].split('\\')
            nameFolder = 'Join_of_{}_files'.format(len(self.dataFolder))
            savepath = '\\'.join(listFolder[0:-2]) + '\\{}\\'.format(nameFolder)
        
        nameExtension = 'R{}_naver{}'.format(self.rudat, self.n_aver)
        
        if saveOnlyProcessed:
            self.saveProcessedData(savepath=savepath+'Processed_data\\', savename='Processed_data_{}.dat'.format(nameExtension))
            return True
        
        print('')
        
        self.saveProcessedData(savepath=savepath+'Processed_data\\', savename='Processed_data_{}.dat'.format(nameExtension))
        self.saveQNDness(save=saveResult, savepath=savepath+self.nameQNDfolder+'\\', savename='QNDness_{}.dat'.format(nameExtension))
    
        print('')
        
        self.plotMeanProba(save=saveResult, savepath=savepath+'Proba_traj\\', savename='Mean_Probabilities_{}.png'.format(nameExtension))
        self.plotMeanTraject(save=saveResult, savepath=savepath+'Stat_traj\\', savename='Statistics_traj_{}.png'.format(nameExtension))
        
        self.plotHistogramBigTraject(save=saveResult, savepath=savepath+'Hist_big_traj\\', savename='Hist_big_traj_{}.png'.format(nameExtension))
        self.plotHistogram_ge(save=saveResult, savepath=savepath+'Hist_big_traj\\', savename='Hist_eg_state_{}.png'.format(nameExtension))
                
        self.plotHistogramDuration(save=saveResult, savepath=savepath+'Hist_dtj\\', savename='Hist_dtj_{}.png'.format(nameExtension))
        self.plotHistogramDuration(log=True, save=saveResult, savepath=savepath+'Hist_dtj\\', savename='Hist_dtj_{}_log.png'.format(nameExtension))
        self.plotCumulativeFrequency(save=saveResult, savepath=savepath+'Hist_dtj\\', savename='CumFreq_dtj_{}.png'.format(nameExtension))
                
        for type_trj in self.listTypeTraject:
            self.colorplotAllTraject(type_trj, save=saveResult, savepath=savepath+'Colorplot\\', savename='Color_all_{}_{}.png'.format(type_trj, nameExtension))
        self.plotAllTraject(percent_plot, savepath=savepath+'Single_traj_naver={}\\'.format(self.n_aver), type_name=nameExtension)
    
        if self.n_aver == 0:
            self.plotDistance_eg(save=saveResult, savepath=savepath, savename='Distance_eg{}.png'.format(nameExtension))
    


