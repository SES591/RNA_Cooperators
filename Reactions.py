######################################################################################################################################
## Copyright Sara Imari Walker 2010
##           Cole Mathis 2015   
######################################################################################################################################
import numpy as np
import random
import math

from Parameters import *
from Polymers import *

######################################################################################################################################
def spontaneous_assembly(fragments, assemblies):
    """Spontaneous Polymer Assembly"""

    import Parameters
    
    Parameters.ks_events += 1
    diceroll = random.random()*Parameters.Ap_s
    checkpoint = 0

    found_ID1 = None
    found_ID2 = None
    for ID1, frag1 in enumerate(fragments):
        
        if frag1.seq == 'W':
            for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'XYZ':
                    # Spontaneous Assembly
                    checkpoint += ks*frag1.tot_count*frag2.tot_count
                    if checkpoint >= diceroll:
                        found_ID2 = ID2
                        found_ID1 = ID1
                        #print found_ID1, found_ID2
                        break
                


        elif frag1.seq == 'WX':
             for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'YZ':
                    # Spontaneous Assembly
                    checkpoint += ks*frag1.tot_count*frag2.tot_count 
                    if checkpoint >= diceroll:
                        found_ID2 = ID2
                        found_ID1 = ID1
                        
                        break
                

        elif frag1.seq == 'WXY':
             for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'Z':
                    # Spontaneous Assembly
                    checkpoint += ks*frag1.tot_count*frag2.tot_count 
                    if checkpoint >= diceroll:
                        found_ID2 = ID2
                        found_ID1 = ID1
                        #print found_ID1, found_ID2
                        break
               
        if found_ID1 != None and found_ID2 != None:
            break
    
    Parameters.fragments, Parameters.assemblies = assemble(found_ID1, found_ID2, fragments, assemblies)
######################################################################################################################################
def assemble(ID1, ID2, fragments, assemblies):
    
    new_igs = fragments[ID1].igs
    new_tag = fragments[ID1].tag
    new_seq = fragments[ID1].seq + fragments[ID2].seq

    found_assembly = False
    for assembly in assemblies:
        if assembly.seq == new_seq and assembly.igs == new_igs and assembly.tag == new_tag:
            assembly_ID = assembly.ID
            found_assembly = True
        if found_assembly == True:
            break

    assemblies[assembly_ID].tot_count += 1.0
    fragments[ID1].tot_count -= 1.0
    fragments[ID2].tot_count -= 1.0
    if fragments[ID1].tot_count < 0.0 or fragments[ID2].tot_count < 0.0:
        print ID1, fragments[ID1].seq, fragments[ID1].tot_count
        print ID2, fragments[ID2].seq, fragments[ID2].tot_count
        exit()
    
    

    return fragments, assemblies

######################################################################################################################################
def catalytic_assembly(fragments, assemblies):
    """Catalyzed Assembly"""

    import Parameters
    diceroll = random.random()*Parameters.Ap_ca
    checkpoint = 0
    #print diceroll
    found_ID1 = None
    found_ID2 = None
    # print [f for f in enumerate(fragments)]
    # raw_input("enter")
    for ID1, frag1 in enumerate(fragments):
    
        if frag1.seq == 'W':
            for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'XYZ':
                    # Catalyzed Assembly
                    for assembly in assemblies:
                        ka = ka_dict[(assembly.igs, frag1.tag)]
                        checkpoint += ka*assembly.tot_count*frag1.tot_count*frag2.tot_count
                        # print "diceroll, checkpoint, assembly count, frag1 %s, fra2 %s" % (frag1.seq, frag2.seq)
                        # print diceroll, checkpoint, assembly.tot_count, frag1.tot_count, frag2.tot_count
                        if checkpoint >= diceroll:
                            # print "Found"
                            found_ID2 = ID2
                            found_ID1 = ID1
                            # raw_input("enter")
                            break
                if found_ID1 != None and found_ID2 != None:
                    break



        elif frag1.seq == 'WX':
             for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'YZ':
                    # Catalyzed Assembly
                    for assembly in assemblies:
                        ka = ka_dict[(assembly.igs, frag1.tag)]
                        checkpoint += ka*assembly.tot_count*frag1.tot_count*frag2.tot_count
                        if checkpoint >= diceroll:
                            found_ID2 = ID2
                            found_ID1 = ID1
                            # raw_input("enter")
                            break
                if found_ID1 != None and found_ID2 != None:
                    break

               

        elif frag1.seq == 'WXY':
             for ID2, frag2 in enumerate(fragments):
                if frag2.seq == 'Z':
                    # Catalyzed Assembly
                    for assembly in assemblies:
                        ka = ka_dict[(assembly.igs, frag1.tag)]
                        checkpoint += ka*assembly.tot_count*frag1.tot_count*frag2.tot_count
                        if checkpoint >= diceroll:
                            found_ID2 = ID2
                            found_ID1 = ID1
                            # raw_input("enter")
                            break
                if found_ID1 != None and found_ID2 != None:
                    break  
        if found_ID1 != None and found_ID2 != None:
            break
    # raw_input("enter")
    Parameters.fragments, Parameters.assemblies = assemble(found_ID1, found_ID2, fragments, assemblies)
               