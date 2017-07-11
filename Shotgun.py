import json
import numpy as np
import pandas as pd
import time
from timeit import default_timer as timer
import sys
import os
from random import randint
import subprocess as sub
from string import Template





class Shotgun(object):
    def __init__(self, molecule, OptType ,functional = ["B3LYP"], basisSet = ["6-311G"]):
        self.OptType = OptType
        #pull out data from state for conveinence?????
        self.molecule = molecule
        #Create Directory for this molecule
        self.directoryName = str(molecule.name)
        if not os.path.exists(self.directoryName):
            os.makedirs(self.directoryName)
        
        # Set up range of basis set and functional
        self.functional = functional
        self.basisSet = basisSet
        
        # Build Index for data frame bases on what options were selected.
        self.size = len(self.functional)*len(self.basisSet)
        self.index = pd.MultiIndex.from_product([self.functional,self.basisSet], names =['Functional','Basis Set'])
        self.results = pd.DataFrame(index = self.index) # Create dataframe of size INDEX
        
        # Fill in data frame with blanks
        self.results['ID'] = np.random.randint(10**5,10**8, size=self.size)
        self.results['Job Status'] =  None
        self.results['Energy'] = np.random.randint(1,1000, size=self.size) # This is just a random energy value
        self.results['Run Time'] = np.nan
        
        if self.size > 10:
            self.maxJobs = 10
        else:
            self.maxJobs = self.size



    def write_input(self, index, row):       
        name = str(row['ID'])
        func = str(index[0])
        bs = str(index[1])
        Optimization  = {'Well': 'TightOpt',
                        'TS': 'TSOpt'}
        
        # Read in input template
        with open('input.template','r') as template:
            temp = template.read()
            s = Template(temp)
        # If template file opens write out input files is a molecular subdirectory
            with open( self.directoryName + '/' + name + '.inp' , 'w') as f:
                out = s.substitute(Name = self.molecule.name,
                                    Level = func,
                                   Basis = bs,
                                   Opt = Optimization[self.OptType],
                                  Charge = self.molecule.charge,
                                  Multi= self.molecule.multiplicity,                                   
                                  XYZ = self.molecule.get_xyz_string())
                f.write(out)
                       
    def write_submit(self, index, row):
        name = str(row['ID'])
        #Write Orca Submit File 
        with open('submit.template' , 'r') as template:
            contents = template.read()
            s = Template(contents)
            substituted = s.substitute( Name = name,
                                        InputName= self.directoryName + '/' + name + '.inp',
                                        OutputName = self.directoryName + '/' + name + '.out') 

            with open('orca.pbs', 'w') as f:
                f.write(substituted)
    
    def submit_job(self,index, row):
        print ("In submit job")
        #write input file
        self.write_input(index,row)
        #write submit file
        self.write_submit(index,row)
        #Submit File
        sub.check_output(['qsub','orca.pbs'])   #submit job
        self.results.set_value(index,'Job Status', 'Running')
     
        
    
    def check_running(self,index, row):
        #Check if a specificly named job is running, techinically returns number of times it is running.
        currentCL = sub.Popen(['qstat | grep -c '+ str(row['ID'])], shell = True, stdout = sub.PIPE)
        running, error = currentCL.communicate()
        print("Check Running get a : " + str(running))
        if int(running) != 0:
            self.results.set_value(index,'Job Status', 'Finished')
        return int(running) #!= 0: string must me cast to int here. Check that at least one instance is queued
 
    def read_output(self):
        pass

    def job_watcher(self):
        jobsToRun = self.size
        runningJobs = 0

        # Main Loop: Run until all jobs are finished.
        while jobsToRun != 0:
            print(self.results)
            time.sleep(30)
            #clean up and tell how many jobs are left to run
            #clear_output()
            sub.Popen(['clear'])
            print("Remaining Jobs: " + str(jobsToRun))
            #Which Molecule are we running calculation for?
            self.molecule.molecule_print()
            
            
            
            #Iterate over DataFrame
            for index, dfRow in self.results.iterrows():
                if dfRow['Job Status'] == 'Running':
                    if self.check_running(index, dfRow): # If Job passes running check
                        print ('RUNNING ::  ' + str(index[0]+" /"+str(index[1])))
                    else:
                        jobsToRun -= 1
                        runningJobs -= 1
                        pass
                
                #Go down dataframe and sumbit empty jobs
                if runningJobs != self.maxJobs:
                    self.submit_job(index, dfRow)
                    print(str(index)+":"+ str(dfRow['Job Status']))
                    runningJobs += 1

    def fire(self, state):        
        print("Firing Calculation")
        print ("-----------------------Entering Submit cycle")
        start = timer()
        self.job_watcher()
        end = timer()
        print ("Total Time : " + str(end-start))
        print ("-----------------------Results")
        print (self.results)
        return self.results
