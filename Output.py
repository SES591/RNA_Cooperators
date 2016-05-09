import math
import numpy as np
######################################################################################################################################
def output_data(exp, infreq, tau,  assemblies, fragments):

    import Parameters
    from Parameters import t_infrequent

    output_fragments(exp, tau, fragments)
    output_assemblies(exp, tau, assemblies)
    coarse_assemblies_rank(exp, tau, assemblies)
    coarse_assemblies_binary(exp, tau, assemblies)


    infreq += t_infrequent
    return infreq
    
######################################################################################################################################
###################### Calculate Ouput ###############################################################################################
######################################################################################################################################
######################################################################################################################################
def output_fragments(exp, t, fragments):
    """Print time-dependent polymer concentrations"""

    import Parameters

    for ID, frag in enumerate(fragments):

        filename = ('%s/%i_fragment_%s.dat' % (Parameters.dirname, exp, ID))
    
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
            
        s = str(t) + '      ' + str(frag.tot_count)
        file.write(s)
        file.write('\n')
        file.close()
######################################################################################################################################
def output_assemblies(exp, t, assemblies):
    """Print time-dependent polymer concentrations for catalytic sequences"""

    import Parameters

    for ID, assembly in enumerate(assemblies):

        filename = ('%s/%i_assembly_%s.dat' % (Parameters.dirname, exp, ID))
    
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
        s = str(t) + '      ' + str(assembly.tot_count)
        file.write(s)
        file.write('\n')
        file.close()
######################################################################################################################################
def coarse_assemblies_binary(exp, t, assemblies):
    """Print time-dependent polymer concentrations for catalytic sequences"""

    import Parameters

    assembly_counts = []
    filename = ('%s/%i_coarse_assemblies.dat' % (Parameters.dirname, exp))
    labels = '0 '
    for ID, assembly in enumerate(assemblies):
        assembly_counts.append(assembly.tot_count)

        if t == 0:
            labels += assembly.igs + assembly.tag + '  '         
    
        
    if(t == 0):
        file = open(filename, 'w')
        file.write(labels)
        file.write('\n')
        normalization = 1.0
    else:
        file = open(filename, 'a') 
        normalization = sum(assembly_counts)
    coarse_fraction = 1.0/len(assembly_counts)
    coarse = []
    for value in assembly_counts:
        mass_fraction = float(value)/normalization
        if mass_fraction >= coarse_fraction:
            coarse.append(1.0)
        else:
            coarse.append(0.0)

    s = str(coarse)
    file.write(s)
    file.write('\n')
    file.close()
######################################################################################################################################
def coarse_assemblies_rank(exp, t, assemblies):
    """Print time-dependent polymer concentrations for catalytic sequences"""

    import Parameters

    assembly_counts = []
    filename = ('%s/%i_rank_coarse_assemblies.csv' % (Parameters.dirname, exp))
    labels = ''
    coarse = ''
    for ID, assembly in enumerate(assemblies):
        assembly_counts.append(assembly.tot_count)

        if t == 0:
            labels += assembly.igs + assembly.tag + ','         
    
        
    if(t == 0):
        label = labels[:-1]
        file = open(filename, 'w')
        file.write(label)
        file.write('\n')
       
    else:
        file = open(filename, 'a') 
        
    sorted_list = sorted(assembly_counts, reverse = True)

    for value in assembly_counts:
        coarse += str(sorted_list.index(value)) + ','

    s = coarse[:-1]
    file.write(s)
    file.write('\n')
    file.close()
######################################################################################################################################
def seq_dictionary(exp, sequences):

    import Parameters
    
    file_name = ('%s/%i_seq_dictionary.dat' % (Parameters.dirname, exp))
    file = open(file_name, 'w')
    

    for ID in range(len(sequences)):
        
        s = str(sequences[ID].seq_ID) + '      ' + str(sequences[ID].sequence)
        file.write(s)
        file.write('\n')

    file.close()
######################################################################################################################################
def diversity(exp, sequences, t):
    """Calculate Shannon Diversity of Extant Population"""

    import Parameters

    I = 0.0
    extant_species = 0.0

    for ID, v in sequences.items():

        if sequences[ID].tot_count != 0:

            extant_species += 1

            p_i = float(sequences[ID].tot_count)/Parameters.Npoly

            I -= p_i*math.log(p_i, 2)


    m_species = float(len(Parameters.monomer_species))

    N = pow(m_species, float(Parameters.R_L))

    max_diversity = -math.log(1.0/N, 2)       #calculates maximum shannon diversity of a population of N = extant_species distinct species       #calculates maximum shannon diversity of a population of N = extant_species distinct species

    if max_diversity != 0:
        I /= max_diversity

    file_name = ('%s/%i_extant_diversity.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
        
    s = str(t) + "        " + str(I)
    file.write(s)
    file.write('\n')

    file.close()

    IG = -I + 1.0   #Calculates information gathered by subtracting off from normalized max diversity -see Krakauer and normalized his measure (0,1) 

    file_name = ('%s/%i_IG.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
        
    s = str(t) + "        " + str(IG)
    file.write(s)
    file.write('\n')

    file.close()

######################################################################################################################################
def rates(exp, t, sequences, monomers):
    """Prints current total propensity for replication, degradation and assembly"""

    import Parameters

    filename = ('%s/%i_rate_Rep.dat' % (Parameters.dirname, exp))
	
    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')
            
    s = str(t) + '      ' + str(Parameters.Ap_r)
    file.write(s)
    file.write('\n')
    file.close()


    filename = ('%s/%i_rate_Death.dat' % (Parameters.dirname, exp))
	
    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')
            
    s = str(t) + '      ' + str(Parameters.Ap_h)
    file.write(s)
    file.write('\n')
    file.close()

    filename = ('%s/%i_rate_SpA.dat' % (Parameters.dirname, exp))
	
    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')
            
    s = str(t) + '      ' + str(Parameters.Ap_s)
    file.write(s)
    file.write('\n')
    file.close()


    if Parameters.k_ch != 0:
        filename = ('%s/%i_rate_catRec.dat' % (Parameters.dirname, exp))
	
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
            
        s = str(t) + '      ' + str(Parameters.Af_ch)
        file.write(s)
        file.write('\n')
        file.close()

    if Parameters.k_cs != 0:
        filename = ('%s/%i_rate_catAssemble.dat' % (Parameters.dirname, exp))
	
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
            
        s = str(t) + '      ' + str(Parameters.Af_cs)
        file.write(s)
        file.write('\n')
        file.close()

    if Parameters.k_ca != 0:
        filename = ('%s/%i_rate_catRecombine.dat' % (Parameters.dirname, exp))
	
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
            
        s = str(t) + '      ' + str(Parameters.Af_ca)
        file.write(s)
        file.write('\n')
        file.close()    

######################################################################################################################################
def output_network_polymers(exp, t, sequences, monomers):
    """Print time-dependent polymer concentrations for catalytic sequences"""


    from Parameters import catalysts, substrates, dirname
    polymers = []
    polymers.extend(catalysts)
    polymers.extend(substrates)
    polymers = list(set(polymers))
    
    for index, ID in enumerate(polymers):

        filename = ('%s/%i_sequence_%s.dat' % (dirname, exp, ID))
    
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
        s = str(t) + '      ' + str(sequences[ID].tot_count)
        file.write(s)
        file.write('\n')
        file.close()
######################################################################################################################################
def output_special_seq(ID, sequences, exp, t):
    """Print concentration of inoculated polymer"""

    import Parameters

    filename = ('%s/%i_sequence_%s.dat' % (Parameters.dirname, exp, sequences[ID].sequence))
	
    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')
    
    s = str(t) + '      ' + str(sequences[ID].tot_count)
    file.write(s)
    file.write('\n')
    file.close()
######################################################################################################################################
def output_all_seqs(sequences, exp, t):
    """Print concentration of all polymers"""

    import Parameters
    for ID, v in sequences.items():
        filename = ('%s/%i_sequence_%s.dat' % (Parameters.dirname, exp, sequences[ID].sequence))
    
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
        
        s = str(t) + '      ' + str(sequences[ID].tot_count)
        file.write(s)
        file.write('\n')
        file.close()
######################################################################################################################################
# Infrequent Calculations
######################################################################################################################################
def population_landscape(exp, t, sequences, monomers):
    """Prints to file the population profile at time t to file"""

    import Parameters

    from Parameters import N_bins

    bin_res = float(Parameters.kr/N_bins)
    
    bin_k = [0]*(N_bins + 1)
    bin_abundance = [0]*(N_bins + 1)


    for i in range(len(bin_k)):

        bin_k[i] = i*bin_res
    
    for ID, v in sequences.items():

        bin_num = int(math.floor((N_bins - 1)*(sequences[ID].k/Parameters.kr)))

        bin_abundance[bin_num] += sequences[ID].tot_count

    fig = plt.figure()   
    ax = fig.add_subplot(111)
    ax.set_ylabel('Abundance')
    ax.set_xlabel('k_value')
    ax.set_xticks(bin_k)
    width = 0.1
   
    rects1 = ax.bar(bin_k, bin_abundance , width, color='cyan', align = 'center')
    
    #plt.setp(plt.gca(), 'yticklabels', [])
    #plt.setp(plt.gca(), 'xticklabels', [])

    file_name = ('%s/%i_population_landscape_%2.f.png' % (Parameters.dirname, exp, t))

    plt.savefig(file_name, format = 'png')    
######################################################################################################################################
def sequence_distribution(exp, t, sequences, monomers):
    """Prints to file the population profile at time t to file"""

    import Parameters

    N_bins = Parameters.R_L + 1
    
    bin_nA = [0]*(N_bins + 1)
    bin_abundance = [0]*(N_bins + 1)

    for i in range(Parameters.R_L + 1):
        """Initialize array to populations"""

        ##############################################
        ### Use binomial coefficient to determine how many sequences have given A/B ratio with: n_A = i, n_B = R_L - i
        ##############################################

        bin_nA[i] = i #number of A monomers in sequence


    for ID, k in sequences.items():

        i = sequences[ID].sequence.count('A')
        bin_abundance[i] += sequences[ID].tot_count
  
 
    fig = plt.figure()   
    ax = fig.add_subplot(111)
    ax.set_ylabel('Abundance')
    ax.set_xlabel('Number of A-monomers in sequence')
    ax.set_xticks(bin_nA)
    width = 0.25
   
    rects1 = ax.bar(bin_nA, bin_abundance , width, color='cyan', align = 'center')
    
    #plt.setp(plt.gca(), 'yticklabels', [])
    #plt.setp(plt.gca(), 'xticklabels', [])

    file_name = ('%s/%i_sequence_distribution_%2.f.png' % (Parameters.dirname, exp, t))

    plt.savefig(file_name, format = 'png')
######################################################################################################################################
def final_data(exp, sequences, mass, runtime):
    """Output final run statistics"""
    import Parameters

    """Write Polymer list to file"""
    filename = ('%s/%i_run_statistics.txt' % (Parameters.dirname, exp))

    
    file = open(filename, 'w')

    file.write("Run parameters \n")
    s = 'Total mass = ' + str(mass) + '\n'
    file.write(s)
    s = 'Length of Sequences = ' + str(Parameters.R_L) + '\n'
    file.write(s)
    s = 'kr = ' + str(Parameters.kr) + '\n'
    file.write(s)
    s = 'kh = ' + str(Parameters.kh) + '\n'
    file.write(s)
    s = 'km = ' + str(Parameters.km) + '\n'
    file.write(s)
    s = 'kc = ' + str(Parameters.kc) + '\n'
    file.write(s)
    s = 'ks = ' + str(Parameters.ks) + '\n'
    file.write(s)
    s = 'mu = ' + str(Parameters.mu) + '\n'
    file.write(s)
    s = 'k_ch = ' + str(Parameters.k_ch) + '\n'
    file.write(s)
    s = 'k_ca = ' + str(Parameters.k_ca) + '\n'
    file.write(s)
    s = 'k_cs = ' + str(Parameters.k_cs) + '\n'
    file.write(s)
    s = 'Run seed = ' + str(Parameters.run_seed) + '\n'
    #file.write(s)
    #s = 'Number of Steps = ' + str(Parameters.tot_step) + '\n' + '\n'
    
    
    file.write("Information about kMC run including run statistics: \n \n")

    file.write("Total run time = ")
    s= str(runtime) + ' sec'
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of replication events = ")
    s = str(Parameters.kr_events)
    file.write(s)
    file.write('\n')

    file.write("Number of error-prone replication events = ")
    s = str(Parameters.mr_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of degradation events =  ")
    s = str(Parameters.kh_events)
    file.write(s)
    file.write('\n')

    file.write("Number of mutagenic events ")
    s = str(Parameters.km_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of cross-over/recombination events ")
    s = str(Parameters.kc_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of spontaneous assembly events ")
    s = str(Parameters.ks_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of catalyzed recycling events ")
    s = str(Parameters.kch_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of catalyzed recombination events ")
    s = str(Parameters.kca_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Number of catalyzed assembly events ")
    s = str(Parameters.kcs_events)
    file.write(s)
    file.write('\n')
    file.write('\n')

    file.write("Final number of monomers = ")
    s = str(Parameters.Nmono)
    file.write(s)
    file.write('\n')

    file.write("Final number of polymers = ")
    s = str(Parameters.Npoly)
    file.write(s)
    file.write('\n')

    file.write("Null Replication events = ")
    s = str(Parameters.null_event)
    file.write(s)
    file.write('\n')

    file.close()

    ''' Save a list of all initial and final replicator composition at the end of the run '''
    initial_composition = [0, 0, 0, 0]
    surviving_composition = [0, 0, 0, 0]

    filename = ('%s/%i_surviving_species.dat' % (Parameters.dirname, exp))
    file = open(filename, 'w')
    
    for ID in range(len(sequences)):
        if sequences[ID].tot_count !=0:

            s = sequences[ID].sequence + '  ' +str(sequences[ID].seq_ID)
            surviving_composition = np.add(sequences[ID].tot_count*np.array(sequences[ID].seq_list), surviving_composition)
            file.write(s)
            file.write('\n')

    file.close()

    for ID in range(len(Parameters.initial_replicators)):
        initial_composition = np.add(Parameters.R_N*np.array(sequences[ID].seq_list), initial_composition)
    compositions = zip(initial_composition, surviving_composition)
    filename = ('%s/%i_compositions.dat' % (Parameters.dirname, exp))
    np.savetxt(filename, compositions)


    if Parameters.output_plots == True:
        import Plotting

        Plotting.output_plots(exp, sequences, Parameters.catalysts, Parameters.substrates)