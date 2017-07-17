from random import randint
from State import State
from Shotgun import Shotgun

import json
import time
import sys
import os

import numpy as np
import pandas as pd
import subprocess as sub
from timeit import default_timer as timer

class Reaction(object):
    def __init__(self, JSON_File = ''):
        self.reaction_states = {}
        ## Reaction reasults as a dictionary with states as KEY, and results dataframe as VALUE
        self.reaction_results= {}
        # Open up the JSON File and quit program if FAIL
        try:
            with open(JSON_File) as file:
                parsed_data = json.load(file)
                self.reactionBase = parsed_data['Base']
        except:
            print("JSON FILE ERROR")
            sys.exit()
        
        self.deltas = parsed_data['Deltas']
        
        #Create States
        for key, state in zip(parsed_data['States'].keys(),parsed_data['States']):
            self.reaction_states[key]= State(parsed_data,key)


        #Create Tree of states _____THIS PART CAN WAIT____

        
    def reaction_print(self):
        for key, state in self.reaction_states.items():
            state.state_print()



    def run_calculations(self, options = None):
        ## TODO Change this to only sumbit each molecule once....
        for key, state in self.reaction_states.items():
            for mol in state.madeUpOf:
                print ("Running Molecule: " + mol.name)
                # maybe just pass default values, but I would maybe refractor
                #  the SHOTGUN CLASS TO EXPECT a LIST
                if options:
                    shotgun = Shotgun(mol, state.Type,options[0],options[1])
                else:
                    shotgun = Shotgun(mol, state.Type)
                # Add to list of results:
                self.reaction_results[mol.name] = shotgun.fire(mol)
                print (self.reaction_results[mol.name])




    def recover_calculations(self, options = None):
        for key, state in self.reaction_states.items():
            for mol in state.madeUpOf:
                if options:
                    shotgun = Shotgun(mol, state.Type,options[0],options[1])
                else:
                    shotgun = Shotgun(mol, state.Type)

                self.reaction_results[mol.name] = shotgun.recover(mol)
                print (self.reaction_results[mol.name])



    def process_results(self):
        #grab an example key
        key = next(iter(self.reaction_results))
        #create blank dataFrame with the correct index
        index  =  self.reaction_results[key].index

        cleanData = pd.DataFrame(index = index)
        for key, state in self.reaction_states.items():
            tempEnergy = pd.DataFrame(index = index)
            tempEnergy[state.key] = 0

            for mol in state.madeUpOf:
                tempEnergy[state.key] =  tempEnergy[state.key] + self.reaction_results[mol.name]["Energy"]
            cleanData[state.key] = tempEnergy[state.key]

        base = pd.DataFrame.copy(cleanData[str(self.reactionBase)])
        for key, state in self.reaction_states.items():
            cleanData[state.key] = cleanData[state.key]-base

        return cleanData


    def subway_graphs(self):
        pass
