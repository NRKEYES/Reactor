from random import randint
from State import State
from Shotgun import Shotgun

import cclib
import json
import time
import sys
import os
import re

import numpy as np
import pandas as pd
import subprocess as sub
from timeit import default_timer as timer

class Reaction(object):
    def __init__(self, JSON_File = ''):
        self.fileName = JSON_File
        self.reaction_states = {}
        ## Reaction reasults as a dictionary with states as KEY, and results dataframe as VALUE
        self.reaction_results= {}
        # Open up the JSON File and quit program if FAIL
        try:
            with open(self.fileName+'.json') as file:
                self.parsedJSON = json.load(file)
                self.reactionBase = self.parsedJSON['Base']
        except:
            print("JSON FILE ERROR")
            sys.exit()
        
        self.deltas = self.parsedJSON['Deltas'] 
        
        #Create States
        for key, state in zip(self.parsedJSON['States'].keys(),self.parsedJSON['States']):
            self.reaction_states[key]= State(self.parsedJSON,key)


        #Create Tree of states _____THIS PART CAN WAIT____

        
    def reaction_print(self):
        for key, state in self.reaction_states.items():
            state.state_print()



    def run_calculations(self, options = ['B3LYP','TVDZ']):
        ## TODO Change this to only sumbit each molecule once....
        for key, state in self.reaction_states.items():
            for mol in state.madeUpOf:
                print ("Running Molecule: " + mol.name)

                shotgun = Shotgun(mol, state.type , directoryName = self.fileName, functional = options[0], basisSet =options[1])

                # Add to list of results:
                self.reaction_results[mol.name] = shotgun.fire(mol)
                print (self.reaction_results[mol.name])




    def recover_calculations(self, options =  ['B3LYP','TVDZ']):
        for key, state in self.reaction_states.items():
            for mol in state.madeUpOf:

                shotgun = Shotgun(mol, state.type , directoryName = self.fileName,functional = options[0], basisSet =options[1])


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



    def shotgun_plots(self):
        #look at the known values and produce shotgun plots comparing error
        # Energy deltas (compare two states)
        # Molecular Dissociations (not states, just molecules)
        # Geometries.
        # Rotational constants (a specific molecule, similar to a  comparison)
        # Frequency
        pass



    def get_output_file_name(self, mol, options =  ['B3LYP','TVDZ']):
        index = pd.MultiIndex.from_product([options[0],options[1]])
        outputName = str(index[0])+str(index[1]+mol.name)
        
        directoryName = self.fileName + '/'+str(mol.name)
        
        outputFile = directoryName + '/' + outputName + '.out'
        return outputFile


    def update_geom(self, options =  ['B3LYP','TVDZ']):
        for key, state in self.reaction_states.items():
            for mol in state.madeUpOf:
                data = self.parsedJSON['Molecules'][mol.name]['xyz']
                print (self.get_output_file_name(mol))
                sys.exit()
                outputFile = cclib.ccopen(self.get_output_file_name(mol))
                try:
                    parsedOutputFile = outputFile.parse()
                    newXYZ = parsedOutputFile.atomcoords[-1]
                except:
                    print ("To picky")
                    sys.exit()
    

                for index, line in enumerate(data):
                    line =  re.sub('\d', '', line) # Clear out numbers
                    line =  re.sub('-', '', line)# Clean out negatives
                    line =  re.sub('\.', '*', line) # swap decimals for *
                    # At  this point just the * places are left like markers for the 

                    line = str.replace(line, '*', str(newXYZ[index][0]), 1)
                    line = str.replace(line, '*', str(newXYZ[index][1]), 2)
                    line = str.replace(line, '*', str(newXYZ[index][2]), 3)

                    data[index] = line

                self.parsedJSON['Molecules'][mol.name]['xyz'] = str(data)

        with open('data.txt', 'w') as outfile:
            json.dump( self.parsedJSON, outfile, indent=4, sort_keys=True)
