from random import randint
from State import State
from Shotgun import Shotgun
from timeit import default_timer as timer


import json
import time
import sys
import os

import numpy as np
import pandas as pd
import subprocess as sub

class Reaction(object):
    def __init__(self, JSON_File = ''):
        self.reaction_states = []
        self.reaction_results= {}
        
        
        # Open up the JSON File and quit program if FAIL
        try:
            with open(JSON_File) as file:
                parsed_data = json.load(file)
        except:
            print("JSON FILE ERROR")
            sys.exit()
        
        #Create States
        for key, state in zip(parsed_data['States'].keys(),parsed_data['States']):
            self.reaction_states.append(State(parsed_data,key))

     
        
        #Create Tree of states _____THIS PART CAN WAIT____

        
    def reaction_print(self):
        for state in self.reaction_states:
            state.state_print()


    def recover_calculation(self):
        pass


    def run_calculations(self, options = None):
        ## Reaction reasults as a dictionary with states as KEY, and results dataframe as VALUE
        recation_results= {}
        
        
        for state in self.reaction_states:
            for mol in state.madeUpOf:
                print ("Running Molecule: " + mol.name)
                # maykbe just pass default values, but I would maybe refractor
                #  the SHOTGUN CLASS TO EXPECT a LIST
                if options:
                    shotgun = Shotgun(mol, state.Type,options[0],options[1])
                else:
                    shotgun = Shotgun(mol, state.Type)
                # Add to list of results:
                self.reaction_results[state.key] = shotgun.fire(mol)
