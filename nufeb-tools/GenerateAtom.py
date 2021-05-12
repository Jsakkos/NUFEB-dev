import random
import argparse
import numpy as np
import pickle
from string import Template
import json # For dealing with metadata
import os # For file level operations
import time # For timing demonstrations
import datetime # To demonstrate conversion between date and time formats
from datafed.CommandLib import API

# arguments to modify the conditions of the simulation seeding
parser = argparse.ArgumentParser(description='Create atom definition files')
parser.add_argument('--n', dest='num', action='store',
                   default=1,
                   help='Create atom definition files for NUFEB with --n #files desired (default is 1)')
parser.add_argument('--r', dest='reps', action='store',
                   default=1,
                   help='Number of replicates')
parser.add_argument('--c',dest='culture_type',action='store',default='co',
                    help='Set culture conditions with --c (co-culture), --ax-c (cyano), --ax-e (e.coli)')
parser.add_argument('--co2', dest='co2', action='store',
                   default=1e3,
                   help='Set initial CO2 concentration (mM)')
parser.add_argument('--d', dest='dims', action='store', type=str,
                   default='1e-4,1e-4,1e-5',
                   help='Set simulation box dimensions (m)')
parser.add_argument('--t', dest='timesteps', action='store',
                   default=35000,
                   help='Number of timesteps to run')
parser.add_argument('--suc', dest='sucrose', action='store',
                   default=1e-19,
                   help='Set initial sucrose concentration (mM)')
parser.add_argument('--grid', dest='grid', action='store',
                   default=2,
                   help='Diffusion grid density (um/grid)')
parser.add_argument('--mono', dest='monolayer', action='store',
                  default=True,
                  help='Set seed generation to monolayer of cells')
parser.add_argument('--u', dest='user', action='store',
                  help='CADES/CNMS user ID')

args = parser.parse_args()

# Simulation box dimensions (m)
Dimensions = [float(x) for x in args.dims.split(',')]

## Species-specific parameters

CellInfo = {'cyano': {'GrowthRate' : round(0.06/3600,7),
       'min_length' : 1e-6, 'max_length' : 5e-6, 'Diameter' : 1e-6, 'Density' : 370,
       'Inertia' : {'ixx' : 0, 'iyy' : 0, 'izz' : 9.2e-23, 'ixy' : 0, 'ixz' : 0, 'iyz' : 0},
         'K_s' : {'light' : 3.5e-4,'o2' : 2e-4, 'suc' : 1e-2,'co2' : 1.38e-4},
        'Yield' : 0.55,'Maintenance' : 0,'Decay' : 0},
         'ecw': {'GrowthRate' : 2.7e-04,
       'min_length' : 1.94e-6, 'max_length' : 2.72e-6, 'Diameter' : 0.73e-6,'Density' : 236,
       'Inertia' : {'ixx' : 0, 'iyy' : 0, 'izz' : 9.2e-23, 'ixy' : 0, 'ixz' : 0, 'iyz' : 0},
         'K_s' : {'light' : 0,'o2' : 1e-3, 'suc' : 3.6,'co2' : 5e-2},
        'Yield' : 0.43,'Maintenance' : 9.50e-7,'Decay' : 2e-5}
}


class Nutrient:
    
    def __init__(self, c, d,mw,state,bc):
        self.concentration = c
        self.diffusion = d
        self.moleculuarWeight = mw
        self.state = state
        self.boundary = bc
        if self.moleculuarWeight is not None:
            self.concentrationNufeb = np.format_float_scientific(self.concentration*self.moleculuarWeight*1e-3,precision=1)
        else:
            self.concentrationNufeb = np.format_float_scientific(self.concentration,precision=1)
        
        
light = Nutrient(1e-1,None,None,'g','nn')
co2 = Nutrient(float(args.co2),1.9e-09,44.01,'l','nn')
o2 = Nutrient(0.28125,2.30e-9,32,'l','nn')
sucrose = Nutrient(float(args.sucrose),5.2e-10,342.3,'l','nn')
gco2 = Nutrient(0,None,44.01,'g','nn')

class Cell:
    
    def __init__(self,Species,Group,Info = CellInfo,dims = Dimensions):
        self.group = Group
        self.species = Species
        self.growth = Info[self.species]['GrowthRate']
        self.min_length = Info[self.species]['min_length']
        self.max_length = Info[self.species]['max_length']
        self.diameter = Info[self.species]['Diameter']
        self.density = Info[self.species]['Density']
        self.Ks = Info[self.species]['K_s']
        self.yld = Info[self.species]['Yield']
        self.maintenance = Info[self.species]['Maintenance']
        self.decay = Info[self.species]['Decay']
        self.inertia = Info[self.species]['Inertia']
        self.boundaries = dims
        self.length = random.uniform(self.min_length,self.max_length)
    def Atom(self):
        self.x = np.format_float_scientific(random.uniform(0,self.boundaries[0]),precision=1)
        self.y = np.format_float_scientific(random.uniform(0,self.boundaries[1]),precision=1)
        self.z = np.format_float_scientific(random.uniform(0,self.boundaries[2]),precision=1)
        return ' '.join(map(str, [self.group, 1, self.density, self.x,self.y,self.z,1]))
    

    def rotate(self,z_dim = False):
        angle = random.uniform(1,180)
        x_displacement = np.format_float_scientific(self.length/2*np.cos(angle*np.pi/360),precision=1)
        y_displacement = np.format_float_scientific(self.length/2*np.sin(angle*np.pi/360),precision=1)
        zd = z_dim
        if zd == True:
            z_angle = random.uniform(1,180)
            z_displacement = np.format_float_scientific(self.length/2*np.sin(z_angle*np.pi/360),precision=1)
        else:
            z_displacement = 0
        return [x_displacement, y_displacement, z_displacement]

    def Bacillus(self):
        return ' '.join(map(str, list(self.inertia.values()) + self.rotate() + [self.diameter]))


# check for runs folder
if not os.path.isdir('runs'):
    os.mkdir('runs')

for n in range(1,int(args.num)+1):
    SucRatio = round(random.random(),3)
    SucPct = int(SucRatio*100)
    if args.culture_type == 'co':
        cell_types = ['cyano','ecw']
        n_cyanos = int(random.uniform(1,100))
        n_ecw = int(random.uniform(1,100))
        n_cells = n_cyanos + n_ecw
        cellCount = {'cyano' : n_cyanos,'ecw' : n_ecw}
        cyGroup = 'group CYANO type 1'
        ecwGroup = 'group ECW type 2'
        cyDiv = f'fix d1 CYANO divide 100 v_EPSdens v_divDia1 {random.randint(1,1e6)}'
        ecwDiv = f'fix d2 ECW divide 100 v_EPSdens v_divDia2 {random.randint(1,1e6)}'
    elif args.culture_type == 'ax-c':
        cell_types = ['cyano']
        n_cyanos = int(random.uniform(1,100))
        n_ecw = 0
        n_cells = n_cyanos
        cellCount = {'cyano' : n_cyanos}
        cyGroup = 'group CYANO type 1'
        ecwGroup = ''
        cyDiv = f'fix d1 CYANO divide 100 v_EPSdens v_divDia1 {random.randint(1,1e6)}'
        ecwDiv = ''
    elif args.culture_type == 'ax-e':
        cell_types = ['ecw']
        n_ecw = int(random.uniform(1,100))
        n_cyanos=0
        n_cells = n_ecw
        cellCount = {'ecw' : n_ecw}

        cyGroup = ''
        ecwGroup = 'group ECW type 1'
        cyDiv = ''
        ecwDiv = f'fix d2 ECW divide 100 v_EPSdens v_divDia2 {random.randint(1,1e6)}'


    NutesNum = len(InitialConditions['Nutrients']['Concentration'])
    for r in range(1,int(args.reps)+1):
        L = [' NUFEB Simulation\r\n\n',f'     {n_cells} atoms \n',
             f'     {len(cell_types)} atom types \n',f'     {NutesNum} nutrients \n\n',
             f'  0.0e-4   {InitialConditions["Dimensions"][0] :.2e}  xlo xhi \n',f'  0.0e-4   {InitialConditions["Dimensions"][1] :.2e}  ylo yhi \n',
             f'  0.0e-4   {InitialConditions["Dimensions"][2] :.2e}  zlo zhi \n\n', ' Atoms \n\n'
             ]

        j = 1
        for c, CellType in enumerate(cell_types,start=1):
            for i in range(j,InitialConditions[CellType]['StartingCells']+j):
                length = random.uniform(InitialConditions[CellType]['min_length'],
                                      InitialConditions[CellType]['max_length'])
                x = random.uniform(0+length,InitialConditions['Dimensions'][0]-length)
                y = random.uniform(0+length,InitialConditions['Dimensions'][1]-length)
                z = random.uniform(0+length,InitialConditions['Dimensions'][2]-length)
                L.append(f'     %d {c} {length :.2e}  {InitialConditions[CellType]["Density"]} {x :.2e} {y :.2e} {z :.2e} {length :.2e} \n'% (i))
                j += 1

        L.append('\n')
        L.append(' Nutrients \n\n')
        for i,nute in enumerate(InitialConditions['Nutrients']['Concentration'].keys()):
            L.append(f'     %d {nute} {InitialConditions["Nutrients"]["State"][nute]} {InitialConditions["Nutrients"]["xbc"][nute]} {InitialConditions["Nutrients"]["ybc"][nute]} {InitialConditions["Nutrients"]["zbc"][nute]} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} \n'% (i+1))

        L.append('\n')
        L.append(' Type Name \n\n')
        for c, CellType in enumerate(cell_types,start=1):
            L.append(f'     {c} {CellType}  \n')
        L.append('\n')
        L.append(' Diffusion Coeffs \n\n')
        for key in InitialConditions['Diff_c'].keys():
            L.append(f'     {key} {InitialConditions["Diff_c"][key]} \n')
        L.append('\n')
        L.append(' Growth Rate \n\n')
        for CellType in cell_types:
            L.append(f'     {CellType} {InitialConditions[CellType]["GrowthRate"]} \n')
        L.append('\n')
        L.append(' Ks \n\n')
        for CellType in cell_types:
            k = f'     {CellType}'
            for key in InitialConditions[CellType]['K_s'].keys():
                k = k + ' ' + str(InitialConditions[CellType]['K_s'][key])
            k = k + f' \n'
            L.append(k)
        L.append('\n')
        for key in InitialConditions["cyano"]['GrowthParams'].keys():
            L.append(' ' + key + f' \n\n')
            for CellType in cell_types:
                L.append(f'     {CellType} {InitialConditions[CellType]["GrowthParams"][key]} \n')
            L.append('\n')


        L.append('\n\n')
        x = int(Dimensions[0]*1e6)
        y = int(Dimensions[1]*1e6)
        z = int(Dimensions[2]*1e6)
        # make atom definition file
        for r in range(1,int(args.reps)+1):
            L = [' NUFEB Simulation\r\n\n',f'     {n_cells} atoms \n',
                 f'     {len(cell_types)} atom types \n',f'     {n_cells} bacilli \n\n',
                 f'  0.0e-4   {x :.2e}  xlo xhi \n',f'  0.0e-4   {y :.2e}  ylo yhi \n',
                 f'  0.0e-4   {z :.2e}  zlo zhi \n\n', ' Atoms \n\n'
                 ]

            # Create list of atoms and bacilli for atom definition file
            atoms_list = []
            bacilli_list = []
            j = 1
            for c, CellType in enumerate(cell_types,start=1):
                cell = Cell(CellType,c)
                for i in range(j,cellCount[CellType]+j):
                    atoms_list.append(str(j) + ' ' + cell.Atom() + f' \n')
                    bacilli_list.append(str(j) + ' ' + cell.Bacillus() + f' \n')
                    j += 1
            atoms = L+atoms_list
            atoms.append('\n')
            atoms.append(' Bacilli \n\n')
            atoms = atoms + bacilli_list
            #write atom definition file
            f= open(f"./runs/atom_{n_cyanos}_{n_ecw}_{SucPct}_{r}.in","w+")
            f.writelines(atoms)


    #write initial conditions pickle file
    dumpfile = open(f"./runs/run_{n_cyanos}_{n_ecw}_{SucPct}.pkl",'wb')
    pickle.dump(InitialConditions,dumpfile)
    dumpfile.close()
    #write Inputscript
    #open the file
    filein = open( './templates/Inputscript.txt' )
    #read it
    src = Template( filein.read() )
    #do the substitution
    result = src.safe_substitute({'n' : n, 'SucRatio' : SucRatio, 'SucPct' : SucPct,
                                  'n_cyanos' : n_cyanos, 'n_ecw' : n_ecw,
                                  'Replicates' : args.reps,'Timesteps' : args.timesteps,
                                  'CYANOGroup' : cyGroup,
                                  'ECWGroup' : ecwGroup,
                                  'Zheight' : InitialConditions["Dimensions"][2],
                                 'CYANODiv'  : cyDiv, 'ECWDiv' : ecwDiv,
                                 'GridMesh' : f'{int(InitialConditions["Dimensions"][0]*1e6/int(args.grid))} {int(InitialConditions["Dimensions"][1]*1e6/int(args.grid))} {int(InitialConditions["Dimensions"][2]*1e6/int(args.grid))}'})
    f= open(f"./runs/Inputscript_{n_cyanos}_{n_ecw}_{SucPct}.lammps","w+")
    f.writelines(result)




    #create DataFed collection to hold the results
    df_api = API()
    df_api.setContext('p/eng107')
    collectionName = f'NUFEB_{n_cyanos}_{n_ecw}_{SucPct}_{x}_{y}_{z}'
    parent_collection = df_api.getAuthUser().split('/')[1]
    coll_msg = df_api.collectionCreate(collectionName,
                                    parent_id=parent_collection)
    global_coll_id = coll_msg[0].coll[0].id

    #write slurm script
    #open the file
    filein = open( './templates/Slurm.txt' )
    #read it
    src = Template( filein.read() )
    #do the substitution
    result = src.safe_substitute({'n' : n, 'job' : f"NUFEB_cyano{n}",
                                    'USER' : args.user,'Replicates'  : args.reps,
                                    'SucPct' : SucPct,'n_cyanos' : n_cyanos,
                                    'n_ecw' : n_ecw,'id': global_coll_id})
    f= open(f"./runs/Inputscript_{n_cyanos}_{n_ecw}_{SucPct}.slurm","w+")
    f.writelines(result)
    #write local run script
    #open the file
    filein = open( './templates/local.txt' )
    #read it
    src = Template( filein.read() )
    #do the substitution
    result = src.safe_substitute({'n' : n, 'SucRatio' : SucRatio, 'SucPct' : SucPct,
                                  'n_cyanos' : n_cyanos, 'n_ecw' : n_ecw,
                                  'Reps' : args.reps,'id': global_coll_id})
    f= open(f"./runs/local_{n_cyanos}_{n_ecw}_{SucPct}.sh","w+")
    f.writelines(result)
