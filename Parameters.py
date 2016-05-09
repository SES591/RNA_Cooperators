#Python Modules
from math import sqrt, pi, pow
import os
import random

#####################################################################################################################################
run_seed = 59   #300
start_run = 0      #start number for experimental number
total_runs = 10  #Total Number of experiments to be run

diverse_runs = 1  #If 1, will generate new random network each exp, else will use the same network each run

"""System Parameters & Output Specifications"""
tau_max = 100.0                                #max simulation time (in dimensionless time units)
t_frequent = 1.0                               #size of time step for early times
t_infrequent = 10.0                             #size of time step for large times
N_bins = 10                                    #Total number of bins for population k_value histogram

output_all = 1              #If = 0 no polymer concentrations output, if = 1 polymer concentrations output for all polymers at each time step
output_plots = True

"""Replicator statistics"""
network = 'Cooperate'               # Network topology

init_frags = 500                    # Initial quantity of each unique sequence in the system
init_assemblies = 1

"""Microscopic Reaction Rates"""
#kh = 1.0                    #degradation rate
ksexp = 4
ks = pow(10.0, -ksexp)                  #spontaneous assembly rate
k_ca = 0.0

"""Monomer Species & Statistics"""
monomer_species = ['W', 'X', 'Y', 'Z']
igs_types = ['A', 'C', 'U', 'G']
tag_types = ['A', 'C', 'U', 'G']

#m = len(M_N)



######################################################################################################################################
### Do not modify parameters below this line!
######################################################################################################################################

"""Global variables not to be changed (program will modify)"""

fragements = []
assembled = []
ka_dict = { ('C', 'G'): 0.0415,
            ('A', 'U'): 0.0319,
            ('U', 'A'): 0.0197,
            ('G', 'C'): 0.0125,
            ('G', 'U'): 0.0091,
            ('A', 'C'): 0.0069,
            ('U', 'G'): 0.0049,
            ('U', 'C'): 0.0038,
            ('U', 'U'): 0.0022,
            ('C', 'A'): 0.0020,
            ('C', 'C'): 0.0016,
            ('G', 'G'): 0.0006,
            ('G', 'A'): 0.0005,
            ('A', 'A'): 0.0004,
            ('C', 'U'): 0.0004,
            ('A', 'G'): 0.0001,
}
selfish_ids = []
cooperative_ids = []
N_assembled = 0   		# Total number of polymers in system (summed over lattice sites)
N_fragments = 0               # Total number of monomers in system
tot_species = 0         # Total number of replicator species that have ever existed in sim



Atot = 0                #Total global propensity for all events
Ap_h = 0                #Total propensity for degradation
Ap_s = 0                #Total propensity for spontaneous assembly
Af_ca = 0               #Total propensity for catalyzed assembly


kh_events = 0            # Tracks number of degredation events
ks_events = 0            # Tracks number of spontaneous assembly events
kca_events = 0          # Tracks number of catalyzed assembly events

repeat_seq = 0
null_event = 0



dirname = 'data/%s_initfrags_%i_ksexp_%i_runseed_%i' % (network, init_frags, ksexp, run_seed)

if not os.path.exists(dirname):
    os.makedirs(dirname)
