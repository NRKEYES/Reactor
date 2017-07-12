from Molecule import Molecule

class State(object):
    def __init__(self, JSON_data, key = ''):     
        self.key = key
        self.madeUpOf = []
        
        self.comesFrom = JSON_data['States'][key]['ComesFrom']
        
        self.goesTo = JSON_data['States'][key]['GoesTo']
        
        self.Type = JSON_data['States'][key]['Type']
        
        for mol in JSON_data['States'][key]['MadeUpOf']:
            self.madeUpOf.append(Molecule(JSON_data['Molecules'][mol],mol))

    def state_print(self):
        print ("--------------------------------------------------------------------")
        print ("Key: " + self.key)
        print ("Comes from:")
        for reagent in self.comesFrom:
            print (reagent)
        print ("Made up of: ")
        for mol in self.madeUpOf:
            mol.molecule_print()
        
        print ("Goes to:")
        for product in self.goesTo:
            print (product) 
            
        print ("Type: " + str(self.Type))

