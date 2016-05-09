# Main Evolution loop for ChemEvoDiff

#import pp
import sys
import time
from Polymers import *
from Initialize import *
from Reactions import *
from Output import *

import numpy as np
import random
import math


# Main 
def main():
    import Parameters
    
    for exp in range(start_run, start_run + total_runs):
    
        start_time = time.clock()

        """Initialize random number generator for each run"""
        Parameters.run_seed = 100*exp + run_seed      
        random.seed(Parameters.run_seed)

        """Clear variables for each new experimental run"""
        Parameters.N_fragments = 0
        Parameters.N_assembled = 0

        Parameters.fragments = []

        # Reset simulation clock, freq_counter keeps track for frequent outputs, infrequent for infrequent outputs
        tau = 0
        freq_counter = 0  
        infreq_counter = 0

        """Initialize monomer and polymer populations"""
        #Parameters.monomers = create_monomers(monomers)
        Parameters.fragments = create_fragments()
        Parameters.assemblies = initalize_assemblies()
        update_propensities( Parameters.fragments,  Parameters.assemblies)

        infreq_counter = output_data(exp, infreq_counter, 0, Parameters.assemblies, Parameters.fragments)

        """Checks so each simulation runs a reasonable length of time"""
        tau_max = Parameters.tau_max #initializes tau_max for each run to the maximal value
        
        while tau < tau_max:
            """Main time evolution loop"""
            """Choose Reaction"""
            
            dice_roll = random.random()*Parameters.Atot
            
            if dice_roll < Parameters.Ap_s:
                spontaneous_assembly(Parameters.fragments, Parameters.assemblies)
                #print "spontaneous_assembly"
            elif dice_roll < (Parameters.Ap_s + Parameters.Ap_ca):
                catalytic_assembly(Parameters.fragments, Parameters.assemblies)
                #print "catalytic_assembly"
    
            update_propensities(Parameters.fragments, Parameters.assemblies)

            #############################################################
            """Check Mass Conservation - If violated evolution loop will exit"""
            # mass = 0
            # mono_count = [0.0]  
            # for ID, seq in Parameters.sequences.items():
                
            #     mass += R_L*seq.tot_count
            #     for m in range(len(mono_count)):
            #         mono_count[m]+= seq.tot_count*seq.seq_list[m]
            # for ID, monomer in Parameters.monomers.items():
            #     if ID == 'W':
            #         m = 0
            #     elif ID =='X':
            #         m = 1
            #     elif ID == 'Y':
            #         m = 2
            #     elif ID == 'Z':
            #         m = 3
            #     mass += monomer.tot_count
            #     mono_count[m] += monomer.tot_count

            # if sum(M_N) != 0 and mass != sum(M_N):
            #     print 'Conservation of mass violated, exiting simulation ... '
            #     break

            # elif sum(M_N) == 0 and mass != R_L*init_sample_size:
            #     print 'Conservation of polymer mass violated, exiting simulation ... '
            #     break
            # elif mono_count != M_N:
            #     print 'Conservation of monomer mass violated, exiting simulation ... '
            #     break
            #############################################################
            """ 3. Adjust Time """

            dice_roll = random.random()
            infreq_counter = output_data(exp, infreq_counter, tau, Parameters.assemblies, Parameters.fragments)
            
            if Parameters.Atot == 0:
                print "Atot == 0"
                break

            else:
                
                tau -= math.log(dice_roll)/Parameters.Atot
                print tau
                #print "\n  \n   \n   \n "
            

            #if tau > freq_counter:
            


        
       
        print 'Saving final data for run number ' + str(exp) + ' ...'
        runtime = time.clock() - start_time
        #final_data(exp, Parameters.sequences, mass, runtime)
        #output_data(exp, infreq_counter, tau, Parameters.assemblies, Parameters.fragments)

# End Function Main ######################################################################################################################################

if __name__=="__main__":
    main()
    
