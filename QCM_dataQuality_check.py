import matplotlib
matplotlib.use('Agg') # non-graphical backend is better here
import numpy as np
import itertools
from matplotlib import pyplot as plt
import pandas as pd
import os 
import argparse

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage[version=3]{mhchem}'

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.style'] = 'normal'
matplotlib.rcParams['font.weight'] = 400
matplotlib.rcParams['font.size'] = 17
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True

 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-i", "--filename", help = "Specify input file", type=str, default='data_2check.xlsx')
 
# Read arguments from command line
args = parser.parse_args()


def data_prep(exp_data):
    """
    Function that prepares QCM data for plotting
    Expects inpit as a list of lists [[overtine_number, dissipation], ...]
    """
    exp_data_pairs = list(itertools.combinations(exp_data, 2))
    X = []
    Y = []
    for k, m in exp_data_pairs:
        overtone1 = k[0]
        overtone2 = m[0]
        freq1 = k[1]
        freq2 = m[1]
        if overtone1 > overtone2:
            Y.append(np.sqrt(overtone2)/np.sqrt(overtone1))
            X.append(freq1/freq2)
        else:
            Y.append(np.sqrt(overtone1)/np.sqrt(overtone2))
            X.append(freq2/freq1)
    return X, Y

def data_clean(Data):
    """
    This function cleans up and reformats the input data from
    Excel files. It expects pandas DataFrame from read_excel input
    data format as in the exemplary data_2check.xlsx file
    """

# clean data to take only valid vals
    overtone = [1, 3, 5, 7, 9, 11, 13, 15]
    o_colName = Data.columns[0]
    d_colName = Data.columns[1]
    f_colName = Data.columns[2]
    
    Data_f_ =  list(Data[f_colName].fillna(0).values)
    Data_d_ =  list(Data[d_colName].fillna(0).values)
    Data_o_ =  list(Data[o_colName].values)
    
    Data_freq = [[Data_o_[i], Data_f_[i]] for i,_ in enumerate(overtone) if Data_f_[i] != 0]
    Data_diss = [[Data_o_[i], Data_d_[i]] for i,_ in enumerate(overtone) if Data_d_[i] != 0]
    
    return Data_freq, Data_diss

#%%


if __name__ == "__main__":
    # Turn interactive plotting off
    plt.ioff()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)
    
    # read in exp data
    filename = args.filename
    # reference_filename = 'reference.xlsx' 
    """
    reference data are based on water - heavy water experiment
    QCMD filename is below
    
    Silica-MQ Water-Deuterium Oxide 8-9-17.qsod
    """
    
    Data = pd.read_excel(filename)
    #Data_ref = pd.read_excel(reference_filename)
    
    """ Data_ref is now included 
        Reference data are for H2O/D2O injection   
    """
    
    data_ref = {'Overtone': [1,
                            3,
                            5,
                            7,
                            9,
                            11,
                            13,
                            15,],
                'Dissipation x1e-6': [np.nan, 
                                        25.31,
                                        20.62,
                                        16.55,
                                        14.98,
                                        13.37,
                                        12.37,
                                        np.nan],
                'Frequency, MHz': [np.nan,
                                    -60.19,
                                    -41.62,
                                    -35.11,
                                    -29.78,
                                    -27.59,
                                    -25.18,
                                    np.nan]}
    
    Data_ref = pd.DataFrame(data_ref)
    
    Data_freq, Data_diss = data_clean(Data)
    frequency_data_reference, dissipation_data_reference = data_clean(Data_ref)
    
    # dissipation check
    """
    D_n = 2beta * sqrt(rho eta)/(n^(1/2) f^*)
    
    D_1/D_2 = (n_2)^(1/2) / (n_1)^(1/2)
    
    matching with the reference means only 
    density viscosity damping of first order
    """
    
    X_dat, Y_dat = data_prep(Data_diss)
    X_ref, Y_ref = data_prep(dissipation_data_reference)
    
    plt.figure()
    plt.plot(np.linspace(0, 1, 20), np.linspace(0, 1, 20), c='gray', linestyle=':', linewidth=1)
    plt.scatter(X_dat, Y_dat, edgecolor='red', facecolor='none', marker='^', label='Your Data')
    plt.scatter(X_ref, Y_ref, edgecolor='blue', facecolor='none', marker='o', label='Reference')
    plt.legend()
    plt.xlabel(r'$\sqrt{n_j/n_i}$')
    plt.ylabel(r'$\Delta D_i/\Delta D_j$')
    plt.title('Dissipation scaling check')
    
    plt.savefig('DissipationScaling_report.png', dpi=300)
    plt.close('all')
    
    # frequency check
    
    """
    f_n/n = 1/C * Gamma - beta*(Delta sqrt(rho eta) ) )/(n^(1/2) 
    
    if no adsorption => Gamma = 0, then
    
    f_1/n / f_2/n = (n_2)^(1/2) / (n_1)^(1/2)
    
    Linear scaling here means no adsorption
    Deviation from linear scaling means adsorption???
    and something else
    """
    X_dat, Y_dat = data_prep(Data_freq)
    X_ref, Y_ref = data_prep(frequency_data_reference)
    
    plt.figure()
    plt.plot(np.linspace(0, 1, 20), np.linspace(0, 1, 20), c='gray', linestyle=':', linewidth=1)
    plt.scatter(X_dat, Y_dat, edgecolor='red', facecolor='none', marker='^', label='Your Data')
    plt.scatter(X_ref, Y_ref, edgecolor='blue', facecolor='none', marker='o', label='Reference')
    plt.legend()
    plt.xlabel(r'$\sqrt{n_j/n_i}$')
    plt.ylabel(r'$\Delta f_i/\Delta f_j$')
    plt.title('Frequency scaling check')
    
    plt.savefig('FrequencyScaling_report.png', dpi=300)
    plt.close('all')