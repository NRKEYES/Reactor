from random import randint
import subprocess as sub
from string import Template
from timeit import default_timer as timer

import time
import glob
import sys
import os
import cclib

import numpy as np
import pandas as pd



class Shotgun(object):
    def __init__(self, molecule, OptType ,functional = ["B3LYP"], basisSet = ["6-311G"]):

        if len(molecule.xyz) > 1:
            self.OptType = OptType
        else:
            self.OptType = 'NoOpt'
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
        self.index = pd.MultiIndex.from_product([self.functional,self.basisSet],
                                                names =['Functional','Basis Set'])
        self.results = pd.DataFrame(index = self.index) # Create dataframe of size INDEX
        
        # Fill in data frame with blanks
        self.results['ID'] = np.random.randint(10**5,10**8, size=self.size)
        self.results['Job Status'] =  None
        self.results['Energy'] = np.nan
        self.results['Run Time'] = np.nan
        
        if self.size > 10:
            self.maxJobs = 10
        else:
            self.maxJobs = self.size



    def filename(self, index, row):
        name = str(index[0])+str(index[1]+self.molecule.name)
        return name



    def write_input(self, index, row):
        name = str(row['ID'])
        func = str(index[0])
        bs = str(index[1])
        Optimization  = {'Well': '! TightOpt',
                            'TS': '! OptTS',
                            'NoOpt': ''}
        ## ADD NO OPT option


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
        try:
            os.remove('orca.pbs')
        except:
            print ("No orca file to delete. Continuing.")

        ID = str(row['ID'])
        outputName = self.filename(index, row)
        #Write Orca Submit File 
        with open('submit.template' , 'r') as template:
            contents = template.read()
            s = Template(contents)
            substituted = s.substitute( ID = ID,
                                        InputName= self.directoryName + '/' + ID + '.inp',
                                        OutputName = self.directoryName + '/' + outputName + '.out')
        with open('orca.pbs', 'w') as writingFile:
            writingFile.write(substituted)
        time.sleep(1)
        return


    def submit_job(self, index, row):
        #write input file
        self.write_input(index,row)
        #write submit file
        self.write_submit(index,row)
        #Submit File
        sub.check_output(['qsub','orca.pbs'])   #submit job
        self.results.set_value(index,'Job Status', 'Running')



    def cleanup(self, index, row):
        # TODO
        # Remove non output files after finishing job
        pass


    def check_running(self, index, row):
        #Check if a specificly named job is running,
        # techinically returns number of times it is running
        currentCL = sub.Popen(['qstat | grep -c '+ str(row['ID'])], shell = True, stdout = sub.PIPE)
        running, error = currentCL.communicate()
        bRunning = bool(int(running))
        #print (bRunning)
        if bRunning: #currently running
            return bRunning
        else:
            self.results.set_value(index,'Job Status', 'Finished')
            return bRunning



    def read_output(self, index= None, row= None, filename = None):
        if filename == None:
            filename = self.directoryName + "/" + self.filename(index, row) + ".out"
            # The else is implicit. Filename will just keep its passed value...
        myOutput = cclib.ccopen(filename)
        try:
            parsed = myOutput.parse()
            if self.OptType == "Well" or self.OptType =='TS':
                if parsed.optdone:
                    self.results.set_value(index,'Energy', parsed.scfenergies[-1])
                else:
                    self.results.set_value(index,'Energy', np.nan)
            else:
                self.results.set_value(index,'Energy', parsed.scfenergies[-1])
        except:
            print ("!!!!!!!!!!!!!!!!!!!!!!!! There was a parsing error.!!!!!!!!!!!!!!!!!!!!!!!!")
            print ("Filename: "+ filename)




    def job_watcher(self):
        jobsToRun = self.size
        runningJobs = 0
        #Which Molecule are we running calculation for?
        self.molecule.molecule_print()
        # Main Loop: Run until all jobs are finished.
        while jobsToRun != 0:
            time.sleep(10)
            print("/" + str(jobsToRun), end = '')
            #Iterate over DataFrame
            for index, row in self.results.iterrows():
                if row['Job Status'] == 'Running':
                    if self.check_running(index, row): # If Job passes running check
                        pass
                        #print ('RUNNING ::  ' + str(index[0]+" /"+str(index[1])))
                    else:
                        self.read_output(index = index, row = row)
                        jobsToRun -= 1
                        runningJobs -= 1

                
                #Go down dataframe and sumbit empty jobs
                if runningJobs != self.maxJobs and row['Job Status'] == None:
                    self.submit_job(index, row)
                    print("\n"+str(index)+":"+ str(row['Job Status']))
                    runningJobs += 1



    def recover(self, state):
        print ("Recovering molecular data... " + self.molecule.name)
        searchTerm = self.directoryName +"/*.out"
        files = glob.glob(searchTerm) # recturns a list of all the .out files in a directory
        for output in files:
            self.read_output(filename = output)
        return self.results



    def fire(self, state):
        print("Firing Calculation")
        print ("---------------------------------------------------------------------------------")
        print ("--------------------------------Entering Submit cycle----------------------------")
        print ("---------------------------------------------------------------------------------")
        start = timer()
        self.job_watcher()
        end = timer()
        print ("Total Time : " + str(end-start))
        return self.results
