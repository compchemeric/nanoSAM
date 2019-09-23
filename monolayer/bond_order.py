import mdtraj as md
import numpy as np
import argparse
import pandas
from collections import deque 


def unit_vector(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)


def trans_frac2(atom_coordinates, number_of_atoms, stepsize):
    trans_frac = []
    
    for i in range(0, number_of_atoms-3, stepsize):
        bond1 = atom_coordinates[:,i+1,:]-atom_coordinates[:,i,:]
        bond2 = atom_coordinates[:,i+2,:]-atom_coordinates[:,i+1,:]
        bond3 = atom_coordinates[:,i+3,:]-atom_coordinates[:,i+2,:]
            
        n1 = np.cross(bond1, bond2)
        n2 = np.cross(bond2, bond3)
        
        m1 = np.cross(n1,n2)
        #m1= unit_vector(m1,1)
        #b_2n = unit_vector(bond2,1)
        
        y =  [n1_*n2_ for n1_,n2_ in zip(m1,bond2)]
        y=np.sum(np.asarray(y),axis=1)
        
        x = [n1_*n2_ for n1_,n2_ in zip(n1,n2)]
        x=np.sum(np.asarray(x),axis=1)
        
        phi = np.rad2deg(np.arctan2(y,x))
        trans=np.any((np.vstack((phi>=150,phi<=-150))),axis=0)
        
        trans_frac.append(sum(trans)/(1.0*len(trans)))
    return trans_frac


def lol(trajectory):
    #Get dataframes sorted
    dataframe = trajectory.topology.subset(trajectory.topology.select('type C and (resname TE8 or resname TE16)')).to_dataframe()[0]
    dataframe['cNumber'] = [int(i[1:]) for i in dataframe.name]
    dataframe = dataframe.sort_values(by=['resSeq', 'cNumber'])
    
    o_dataframe = trajectory.topology.subset(trajectory.topology.select('type O and (resname TE8 or resname TE16)')).to_dataframe()[0]
    o_dataframe['oNumber'] = [int(i[1:]) for i in o_dataframe.name]
    o_dataframe = o_dataframe.sort_values(by=['resSeq', 'oNumber'])
    
    #Prepare dataframes for Anchor, OEG and Matrix
    anchor_df = dataframe.loc[dataframe['cNumber'] >= 31]
    oeg_df = dataframe.loc[dataframe['cNumber'] >= 17]
    oeg_df = oeg_df.loc[oeg_df['cNumber'] <= 30]
    matrix_df = dataframe.loc[dataframe['cNumber'] <= 16]
    
    matrix_number_of_c_atoms = 16
    oeg_number_of_c_atoms = 14
    oeg_number_of_o_atoms = 8
    anchor_number_of_c_atoms = 8
    
    matrix_indices = matrix_df.values[:,0].astype(int)
    anchor_indices = anchor_df.values[:,0].astype(int)
    oeg_c_indices = deque(oeg_df.values[:,0].astype(int))
    oeg_o_indices = deque(o_dataframe.values[:,0].astype(int))
    
    
    #Append oeg in order O-C-C
    oeg_indices = []
    
    while oeg_o_indices:
        for _ in range(7):
            oeg_indices.append(oeg_o_indices.popleft())
            oeg_indices.append(oeg_c_indices.popleft())
            oeg_indices.append(oeg_c_indices.popleft())
        
        oeg_indices.append(oeg_o_indices.popleft())
    
    #Get actual coordinates
    matrix_atom_coordinates = trajectory.xyz[:, matrix_indices].reshape((-1, matrix_number_of_c_atoms, 3))
    oeg_atom_coordinates = trajectory.xyz[:, oeg_indices].reshape((-1, oeg_number_of_c_atoms+oeg_number_of_o_atoms,3))
    anchor_atom_coordinates = trajectory.xyz[:, anchor_indices].reshape((-1, anchor_number_of_c_atoms, 3))
    
    anchor_trans_frac = trans_frac2(anchor_atom_coordinates, anchor_number_of_c_atoms, 1)
    oeg_trans_frac = trans_frac2(oeg_atom_coordinates, oeg_number_of_c_atoms+oeg_number_of_o_atoms, 3)
    matrix_trans_frac = trans_frac2(matrix_atom_coordinates, matrix_number_of_c_atoms, 1)
    
    print(np.mean(matrix_trans_frac))
    print(np.mean(oeg_trans_frac))
    print(np.mean(anchor_trans_frac))
          
          
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Lattice constant calculation')
    parser.add_argument('--trajectory', help='Trajectory input path', required=True)
    parser.add_argument('--topology', help='Topology input path', required=True)
    #parser.add_argument('--output', help='Path to output', required=True)
    
    #parser.add_argument('--txt', dest='format', help='Output as txt', action='store_const', const='txt')
    #parser.add_argument('--numpy', dest='format', help='Output as numpy readable file', action='store_const', const='numpy')
    #parser.set_defaults(format='txt')
    
    #parser.add_argument('--stream', dest='stream', help='Enable iterative file stream for input trajectory', action='store_true')
    #parser.add_argument('--no-stream', dest='stream', help='Disable iterative file stream for input trajectory', action='store_false')
    #parser.set_defaults(stream=False)
    
    #parser.add_argument('--chunks', help='Size of chunks to be streamed if streaming is enabled', type=int)
    parser.add_argument('--verbose', dest='verbose', help='Verbose', action='store_const', const=True)
    parser.set_defaults(verbose=False)
    
    parser.add_argument('--frameskip', dest='frameskip', help='Only evaluate every nth frame', type=int)
    parser.set_defaults(frameskip=1)
    
    args = parser.parse_args()
    
    if args.verbose:
        print("Loading trajectory")
        
    trajectory = md.load(args.trajectory,top=args.topology, stride=args.frameskip)
    
    if args.verbose:
        print("Trajectory loaded")
        
    time = trajectory.time
    no_frames = trajectory.n_frames
    
    lol(trajectory)