import json
import numpy as np
import pandas as pd

class Molecule(object):
    def __init__(self, JSON_data, name = ''):
        
        self.name = name
        self.charge = JSON_data['Charge']
        self.multiplicity = JSON_data['Multiplicity']
        
        self.xyz = JSON_data['xyz']  # This is just lines of text.
        
    def get_xyz_string(self):
        xyzString = ""
        for atom in self.xyz:
            xyzString = xyzString + atom + "\n"
        return xyzString
        
    def molecule_print(self):
        print ("    Name: " + self.name)
        print ("    Charge: " + str(self.charge))
        print ("    Multiplicity: " + str(self.multiplicity))
                
        for atom in self.xyz:
            print (u"\u269B  " +atom)
