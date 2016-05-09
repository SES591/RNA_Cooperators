
import numpy as np
######################################################################################################################################
class Assembly(object):
    """A Catalyst Sequence"""
    
    def __init__(self, tot_count, igs, tag, fragment, ID):
        #print("A new polymer has been born!")
        self.tot_count = tot_count
        self.igs = igs
        self.tag = tag
        self.ID = ID
        self.seq = fragment

        self.Ap_h = 0.0
        self.Ap_a = 0.0
        self.Ap_fa = 0.0
        
    def talk(self):
        print("Hi. I'm a polymer.", "\n")
######################################################################################################################################
class Fragment(object):
    """A Fragment Sequence"""
    
    def __init__(self, tot_count, igs, tag, fragment, fragment_ID):
        #print("A new polymer has been born!")
        self.tot_count = tot_count
        self.igs = igs
        self.tag =  tag
        self.ID = fragment_ID
        self.seq = fragment

        self.Ap_h = 0.0
        self.Ap_a = 0.0
        self.Ap_fa = 0.0
######################################################################################################################################        
class Monomer(object):
    """A Monomer Species"""
    
    def __init__(self, tot_count, species, igs, tag, ID):
        #print("A new polymer has been born!")
        self.tot_count = tot_count
        self.species = species
        self.tag = tag
        self.igs = igs
        self.ID = ID