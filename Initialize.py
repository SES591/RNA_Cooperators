import math
import numpy as np
import random

from Polymers import *
import itertools

def perm_unique(elements):
    """Generates all unique permutations of list in elements"""
    eset=set(elements)
    listunique = [[i,elements.count(i)] for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)
def perm_unique_helper(cl,lst,d):
    if d < 0:
        yield tuple(lst)
    else:
        for (ind,(j,i)) in enumerate(cl):
            if i>0:
                lst[d]=j
                cl[ind][1]-=1
                for g in  perm_unique_helper(cl,lst,d-1):
                    yield g
                cl[ind][1]+=1
###############################################################################
def create_fragments():
    """This polymer list which usually consists only of monomers to start"""

    import Parameters
    import Polymers
    import Output
    from Parameters import ka_dict, init_frags
    from itertools import permutations

    ### Generates all possible permutations with repetition (all possible sequences), created as a list of tuples
    if Parameters.network == 'Cooperate':
        all_frags = []
        ID = 0 #Initializes first ID number

        '''S1 '''
        frag = 'W'
        igs ='U'
        tag = 'A'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID += 1
        Parameters.N_fragments += init_frags

        frag = 'XYZ'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID += 1
        Parameters.N_fragments += init_frags

        '''S2 '''
        frag = 'WX'
        igs = 'A'
        tag = 'U'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        frag = 'YZ'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        '''S3 '''
        frag = 'WXY'
        igs = 'C'
        tag = 'G'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        frag = 'Z'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        '''E1 '''
        frag = 'W'
        igs ='U'
        tag = 'G'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        # frag = 'XYZ'
        # igs = None
        # tag = None
        # all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # # print all_frags[ID].seq,all_frags[ID].tot_count 
        # ID +=1
        # Parameters.N_fragments += init_frags

        '''E2 '''
        frag = 'WX'
        igs = 'A'
        tag = 'A'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        # frag = 'YZ'
        # igs = None
        # tag = None
        # all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # # print all_frags[ID].seq,all_frags[ID].tot_count 
        # ID +=1
        # Parameters.N_fragments += init_frags

        '''E3 '''
        frag = 'WXY'
        igs = 'C'
        tag = 'U'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags



    elif Parameters.network == 'Selfish':
        all_frags = []
        ID = 0 #Initializes first ID number

        '''S1 '''
        frag = 'W'
        igs ='U'
        tag = 'A'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID += 1
        Parameters.N_fragments += init_frags

        frag = 'XYZ'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID += 1
        Parameters.N_fragments += init_frags

        '''S2 '''
        frag = 'WX'
        igs = 'A'
        tag = 'U'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        frag = 'YZ'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        '''S3 '''
        frag = 'WXY'
        igs = 'C'
        tag = 'G'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        frag = 'Z'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

    elif Parameters.network == 'COnly':
        all_frags = []
        ID = 0 #Initializes first ID number

        frag = 'Z'
        igs = None
        tag = None
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        '''E1 '''
        frag = 'W'
        igs ='U'
        tag = 'G'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        # frag = 'XYZ'
        # igs = None
        # tag = None
        # all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # # print all_frags[ID].seq,all_frags[ID].tot_count 
        # ID +=1
        # Parameters.N_fragments += init_frags

        '''E2 '''
        frag = 'WX'
        igs = 'A'
        tag = 'A'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        # frag = 'YZ'
        # igs = None
        # tag = None
        # all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # # print all_frags[ID].seq,all_frags[ID].tot_count 
        # ID +=1
        # Parameters.N_fragments += init_frags

        '''E3 '''
        frag = 'WXY'
        igs = 'C'
        tag = 'U'
        all_frags.append(Fragment(init_frags, igs, tag, frag, ID))
        # print all_frags[ID].seq,all_frags[ID].tot_count 
        ID +=1
        Parameters.N_fragments += init_frags

        

    return all_frags
#####################################################################################################
def initalize_assemblies():
    import Parameters
    import Polymers
    import Output
    from Parameters import ka_dict, init_frags, init_assemblies
    from itertools import permutations

    all_assemblies = []
    selfish_ids = []
    cooperative_ids = []
    ID = 0
    if Parameters.network == 'Cooperate':
        '''S1 '''
        frag = 'WXYZ'
        igs ='U'
        tag = 'A'
        all_assemblies.append( Assembly(0, igs, tag, frag, ID) )
        selfish_ids.append(ID)
        ID += 1


        '''S2 '''
        frag = 'WXYZ'
        igs = 'A'
        tag = 'U'
        all_assemblies.append(Assembly(0, igs,tag, frag, ID))
        selfish_ids.append(ID)
        ID +=1

        '''S3 '''
        frag = 'WXYZ'
        igs = 'C'
        tag = 'G'
        all_assemblies.append(Assembly(0, igs, tag,  frag, ID))
        selfish_ids.append(ID)
        ID +=1

        '''E1 '''
        frag = 'WXYZ'
        igs ='U'
        tag = 'G'
        all_assemblies.append(Assembly(0, igs, tag, frag, ID))
        cooperative_ids.append(ID)
        ID +=1

        '''E2 '''
        frag = 'WXYZ'
        igs = 'A'
        tag = 'A'
        all_assemblies.append(Assembly(0, igs, tag,  frag, ID))
        cooperative_ids.append(ID)
        ID +=1

        '''E3 '''
        frag = 'WXYZ'
        igs = 'C'
        tag = 'U'
        all_assemblies.append(Assembly(0, igs, tag, frag, ID))
        cooperative_ids.append(ID)
        ID +=1

    elif Parameters.network == 'Selfish':
        '''S1 '''
        frag = 'WXYZ'
        igs ='U'
        tag = 'A'
        all_assemblies.append( Assembly(0, igs, tag, frag, ID) )
        selfish_ids.append(ID)
        ID += 1


        '''S2 '''
        frag = 'WXYZ'
        igs = 'A'
        tag = 'U'
        all_assemblies.append(Assembly(0, igs,tag, frag, ID))
        selfish_ids.append(ID)
        ID +=1

        '''S3 '''
        frag = 'WXYZ'
        igs = 'C'
        tag = 'G'
        all_assemblies.append(Assembly(0, igs, tag,  frag, ID))
        selfish_ids.append(ID)
        ID +=1   
        
    elif Parameters.network == 'COnly':
        '''E1 '''
        frag = 'WXYZ'
        igs ='U'
        tag = 'G'
        all_assemblies.append(Assembly(0, igs, tag, frag, ID))
        cooperative_ids.append(ID)
        ID +=1

        '''E2 '''
        frag = 'WXYZ'
        igs = 'A'
        tag = 'A'
        all_assemblies.append(Assembly(0, igs, tag,  frag, ID))
        cooperative_ids.append(ID)
        ID +=1

        '''E3 '''
        frag = 'WXYZ'
        igs = 'C'
        tag = 'U'
        all_assemblies.append(Assembly(0, igs, tag, frag, ID))
        cooperative_ids.append(ID)
        ID +=1

    Parameters.selfish_ids = selfish_ids
    Parameters.cooperative_ids = cooperative_ids


    return all_assemblies


#####################################################################################################
###########    Initialize and Update Reaction Propensities    #######################################
#####################################################################################################
def update_propensities(fragments, assemblies):
    """Initializes stochastic reaction propensities at the start of each day"""

    import Parameters
    from Parameters import ka_dict, ks
    #Initialize local variables
    Parameters.Atot = 0
    Parameters.Af_ca = 0

    Ap_ca = 0.0
    Ap_s = 0.0
    
    for ID1, frag1 in enumerate(fragments):
        
        if frag1.seq == 'W':
            # print frag1.seq
            for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'XYZ':
                    # Spontaneous Assembly
                    Ap_s += ks*frag1.tot_count*frag2.tot_count
                    # print ks, frag1.tot_count, frag2.tot_count, ks*frag1.tot_count*frag2.tot_count
                    # Catalyzed Assembly
                    for assembly in assemblies:
                        ka = ka_dict[(assembly.igs, frag1.tag)]
                        Ap_ca+= ka*assembly.tot_count*frag1.tot_count*frag2.tot_count
                


        elif frag1.seq == 'WX':
             for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'YZ':
                    Ap_s += ks*frag1.tot_count*frag2.tot_count 
                    # Catalyzed Assembly
                    for assembly in assemblies:
                        ka = ka_dict[(assembly.igs, frag1.tag)]
                        Ap_ca+= ka*assembly.tot_count*frag1.tot_count*frag2.tot_count

                

        elif frag1.seq == 'WXY':
            
            for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'Z':
                    Ap_s += ks*frag1.tot_count*frag2.tot_count 
                    # Catalyzed Assembly
                    for assembly in assemblies:
                        ka = ka_dict[(assembly.igs, frag1.tag)]
                        Ap_ca+= ka*assembly.tot_count*frag1.tot_count*frag2.tot_count
                
    
    Parameters.Ap_ca = Ap_ca
    Parameters.Ap_s = Ap_s

    if Parameters.Ap_s < 0:

        print 'Spontaneous Assembly propensity has gone negative'
        exit()

    elif Parameters.Ap_ca < 0:

        print 'Catalyzed assembly propensity has gone negative'
        exit()

    Parameters.Atot = Parameters.Ap_s + Parameters.Ap_ca
    Parameters.fragments = fragments
    Parameters.assemblies = assemblies