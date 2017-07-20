import numpy as np
import pandas as pd

import cclib


h = 6.62608E-34
h2ok = 3.18002e-44 # h^2/k
clight = 2.997925e8 
Tc = 298.15 
RJ = 8.31451 
NA = 6.02214E+23 
pi = 3.141592654 
kB = 1.38066E-23 
Po = 1.0e5 # use Po = 101325 for 1 atm standard state
lamdainv = np.power(2*pi*kB/1000/NA,1.5)/np.power(h,3)*RJ/Po/NA
comment = "" 
group = "forced" 
sigma = 1 # rotational symmetry number=order of rotational subgroup
symtol = 0.02   # tolerance for symmetry matching
symrange = 25.0 # factor for increasing symtol for large molecules
TK = 0 # temperature



h = 6.62608E-34
h2ok = 3.18002e-44 # h^2/k
clight = 2.997925e8 
Tc = 298.15 
RJ = 8.31451 
NA = 6.02214E+23 
pi = 3.141592654 
kB = 1.38066E-23 
Po = 1.0e5 # use Po = 101325 for 1 atm standard state
lamdainv = np.power(2*pi*kB/1000/NA,1.5)/np.power(h,3)*RJ/Po/NA
comment = "" 
group = "forced" 
sigma = 1 # rotational symmetry number=order of rotational subgroup
symtol = 0.02   # tolerance for symmetry matching
symrange = 25.0 # factor for increasing symtol for large molecules
TK = 0 # temperature



def get_mass(element_number):
    mass = (1.00794,4.002602,6.941,9.012182,10.811,12.0107,14.0067,15.9994,
            18.9984032,20.1797,22.98976928,24.305,26.9815386,28.0855,30.973762,
            32.065,35.453,39.948,39.0983,40.078,44.955912,47.867,50.9415,51.9961,
            54.938045,55.845,58.933195,58.6934,63.546,65.39,69.723,72.64,74.9216,78.96,
            79.904,83.798,85.4678,87.62,88.90585,91.224,92.906, 38,95.94,97.9072,101.07,
            102.905, 50,106.42,107.8682,112.411,114.818,118.71,121.76,127.6,126.904, 47,
            131.293,132.9054519,137.327,138.90547,140.116,140.90765,144.242,144.9127,
            150.36,151.964,157.25,158.92535,162.5,164.930, 32,167.259,168.93421,173.04,
            174.967,178.49,180.94788,183.84,186.207,190.23,192.217,195.084,196.966569,
            200.59,204.3833,207.2,208.9804,208.9824,209.9871,222.0176,223.0197,226.0254,
            227.0277,232.03806,231.03588,238.02891,237.0482,244.0642,243.0614,247.0704,
            247.0703,251.0796,252.0830,257.0951,258.0984,259.1010,262.1097,261.1088,
            262,266,264,277,268,271,272,285,284,289,288,292)
    
    return mass[element_number-1]


def moment(natm,a,x,y,z,m):
    # calculate the moment of inertia matrix
    sxx = 0.0 
    syy = 0.0 
    szz = 0.0 
    sxy = 0.0 
    sxz = 0.0 
    syz = 0.0 
    xcm = 0.0 
    ycm = 0.0 
    zcm = 0.0 
    mt = 0.0 
    
    for i in range(0, natm,1):
        xcm += m[i]*x[i]
        ycm += m[i]*y[i]
        zcm += m[i]*z[i]
        mt += m[i] 

    xcm = xcm/mt
    ycm = ycm/mt
    zcm = zcm/mt
    
    for i in range(0, natm,1):
        x[i] = x[i]-xcm
        y[i] = y[i]-ycm
        z[i] = z[i]-zcm

    for i in range(0, natm,1):
        sxx += m[i]*x[i]*x[i]
        syy += m[i]*y[i]*y[i]
        szz += m[i]*z[i]*z[i]
        sxy += m[i]*x[i]*y[i]
        sxz += m[i]*x[i]*z[i]
        syz += m[i]*y[i]*z[i]
        
    a[0][0] = syy+szz
    a[1][0] = -sxy
    a[1][1] = sxx+szz
    a[2][0] = -sxz
    a[2][1] = -syz
    a[2][2] = sxx + syy
    return a
    


    
def get_Rotational_Constant(form,Inert):
    # calculate RC in different units
    str = ""
    if ( Inert>0.001 ):
        Inert = Inert/NA/1000*1e-20 ;
        RC = h/8/pi/pi/Inert/clight/100 ; # cm-1
    else:
        RC = 1000000
    return RC



    

def rotation(file):
    myTest = cclib.ccopen(file)
    try:
        parsed = myTest.parse()
    except:
        print ("To picky")
    
    coords = parsed.atomcoords[-1]




    natm = parsed.natom # number atoms

    Inertia = {'x':0.0,
                'y':0.0,
                'z':0.0}


    # matrix to diagonalize starting at i=1
    # There may be an indexing issue here
    a = np.zeros((3,3))

    MM = 0.0
    # a is based on these bad boys
    x = np.zeros(natm)
    y = np.zeros(natm)
    z = np.zeros(natm)
    m = np.zeros(natm)
    
    
    for i in range(0, natm,1):
        x[i] = coords[i][0]
        y[i] = coords[i][1]
        z[i] = coords[i][2]

        try:
            m[i] = get_mass(parsed.atomnos[i])
            MM += m[i]
        except:
            print("Bad Element")


    a = moment(natm,a, x, y ,z,m)
    np.linalg.eig(a)
    eigenValues , matrix= np.linalg.eigh(a)

    #   Inertia[i] = a[i+(i*i-i)/2] 
    Inertia['x'] = eigenValues[2]
    Inertia['y'] = eigenValues[1]
    Inertia['z'] = eigenValues[0]
    
    result = pd.Series([0.0,0.0,0.0])
    
    # These are the values I actually want.
    result[0] = get_Rotational_Constant(1,Inertia['x'])
    result[1] = get_Rotational_Constant(1,Inertia['y']) 
    result[2] = get_Rotational_Constant(1,Inertia['z'])

    return result

test = rotation('MnD/6MnOH+/M066-31++G**6MnOH+.out')
test
