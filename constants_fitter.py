import numpy as np
import scipy
print(scipy.__version__)
from scipy.optimize import curve_fit

np.random.seed(0)

# get number of available 4D bins 
N=6
n_eta_bins = 2

# 4D bin labels
labels = ['0_0_1_0', '0_0_0_0', '0_0_0_0', '0_0_0_0', '0_0_0_0','0_0_0_0']
index_A_plus = [0 , 0 , 0 , 0, 0, 0]
index_A_minus = [1 , 0 , 0 , 0, 0, 0]
index_e_plus = [x + n_eta_bins for x in index_A_plus]
index_e_minus = [x + n_eta_bins for x in index_A_minus]
index_M_plus = [x + 2*n_eta_bins for x in index_A_plus]
index_M_minus = [x + 2*n_eta_bins for x in index_A_minus]

# for debugging define 2 arrays size N with dummy k values 
k_plus = np.linspace(0.1,1.1,6)
k_minus = np.linspace(1.,2., 6)

# scale^2 data dummy
scale_sq = np.linspace(1.,2., 6)

# params is a 1D array of size 3*number_eta_bins, idices from 0 to n-1 contain A, from n to 2n-1 epsilon, from 2n to 3n-1 M
# Initial parameter guess
params0 = np.zeros(3*n_eta_bins)

# function from label name to give eta,pt,eta,pt bin index
def get_index(label):
    eta_plus_index = int(label.split("_")[0])
    pt_plus_index = int(label.split("_")[1])
    eta_minus_index = int(label.split("_")[2])
    pt_minus_index = int(label.split("_")[3])
    return [eta_plus_index, pt_plus_index, eta_minus_index, pt_minus_index]

# function from pt bin index to give k
def get_k(pt_bin_index):
    # pt binning must match the one used to make the histos
    pt_binning = [25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0]
    # TODO could add the pT histo per eta,pt,eta,pt to get more accurate average pT
    return 2.0 / (pt_binning[pt_bin_index] + pt_binning[pt_bin_index+1]) # k = 1 / avg_pT

# define arrays size N with values of k
#for i in range(N):
#    k_plus[i] = get_k(get_index(labels[i])[1])
#    k_minus[i] = get_k(get_index(labels[i])[3])

def func(index_of_4D_bin, *params): 
    
    # indices to get A,e,M
    index_eta_plus = get_index(labels[index_of_4D_bin])[0] 
    index_eta_minus = get_index(labels[index_of_4D_bin])[2]

    index_A_plus = index_eta_plus
    index_e_plus = index_eta_plus + n_eta_bins
    index_M_plus = index_eta_plus + 2*n_eta_bins
    index_A_minus = index_eta_minus
    index_e_minus = index_eta_minus + n_eta_bins
    index_M_minus = index_eta_minus + 2*n_eta_bins

    kp = 1.0
    km = 1.0 

    # Define v = (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k)
    v = ( 1 + params[index_A_plus] - params[index_e_plus]*kp + params[index_M_plus]/kp )*( 1 + params[index_A_minus] - params[index_e_minus]*km - params[index_M_minus]/km )
    
    return v

index_array = np.array([0,1,2,3,4,5])
#index_array = 0

# Compute solution providing initial guess params0, k input, and scale^2 input
sol, cov = curve_fit(func, index_array, scale_sq, params0)
# TODO match sol to the right parameters
print sol, cov



