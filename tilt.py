import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA
import argparse


#Helper function to split up an array into multiple arrays of fixed length
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]
     

def unit_vector(vector):
  return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
  v1_u = unit_vector(v1)
  v2_u = unit_vector(v2)
  return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def angles(trajectory, resolution):
    selection = []
    
    #Hanlde resolution selection
    if resolution.upper() == "CG":
        selection = trajectory.topology.select('(resname == FIL or resname == TE8) and (name == C1B or name == C2B or name == C3B or name == ACB)')
    elif resolution.upper() == "AA":
        selection = trajectory.topology.select('(resname == FIL or resname == TE8) and (name == C1 or name == C2 or name == C3 or name == C4 or name == C5 or name == C6 or name == C7 or name == C8 or name == C9 or name == C10 or name == C11 or name == C12 or name == C13 or name == C14 or name == C15 or name == C16)')
    else:
        raise Exception('Resolution not valid')
        return
        
    molecule_list = list(chunks(selection, 4))
    
    z_axis = [0,0,1]
    pca = PCA(n_components=1)
    
    avg_molecule_angles = []
    std_devs = []
    
    for frame in range(trajectory.n_frames):
        molecule_angles = []

        for molecule in molecule_list:
            atom_coordinates = trajectory.xyz[frame,molecule]
            pca.fit(atom_coordinates)
            
            vec=pca.components_[0,:]
            vec=vec*np.sign(vec[2])
            vec=vec/(np.sqrt(np.sum(vec**2)))
            
            molecule_angles.append(angle_between(vec, z_axis))
        
        avg_molecule_angles.append(np.mean(molecule_angles))
        std_devs.append(np.std(molecule_angles))
        
    return avg_molecule_angles, std_devs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tilt angle calculation')
    parser.add_argument('--trajectory', help='Trajectory input path', required=True)
    parser.add_argument('--topology', help='Topology input path', required=True)
    parser.add_argument('--output', help='Path to output', required=True)
    
    parser.add_argument('--txt', dest='format', help='Output as txt', action='store_const', const='txt')
    parser.add_argument('--numpy', dest='format', help='Output as numpy readable file', action='store_const', const='numpy')
    parser.set_defaults(format='txt')
    
    parser.add_argument('--stream', dest='stream', help='Enable iterative file stream for input trajectory', action='store_true')
    parser.add_argument('--no-stream', dest='stream', help='Disable iterative file stream for input trajectory', action='store_false')
    parser.set_defaults(stream=False)
    
    parser.add_argument('--chunks', help='Size of chunks to be streamed if streaming is enabled', type=int)
    parser.add_argument('--frameskip', dest='frameskip', help='Only evaluate every nth frame', type=int)
    parser.set_defaults(frameskip=1)
    
    parser.add_argument('--resolution', dest='resolution', help='Coarse grained(CG) or AA')
    parser.set_defaults(resolution='cg')
    
    args = parser.parse_args()
    avg_molecule_angles = []
    std_devs = []
    no_frames = 0
    time_list = []
    
    #Handle streaming or completely loading the file
    if args.stream:
        for chunk in md.iterload(args.trajectory,top=args.topology, chunk = args.chunks, stride=args.frameskip):
            avg_molecule_angles_chunk, std_devs_chunk = angles(chunk, args.resolution)
            avg_molecule_angles.extend(avg_molecule_angles_chunk)
            std_devs.extend(std_devs_chunk)
            time_list.extend(chunk.time)
            no_frames += chunk.n_frames
        
    else:
        trajectory = md.load(args.trajectory,top=args.topology, stride=args.frameskip)
        
        time_list = trajectory.time
        no_frames = trajectory.n_frames
        
        avg_molecule_angles, std_devs = angles(trajectory, args.resolution)
    
    #Save the results
    if args.format != 'numpy':
        with open(args.output, "w") as text_file:
            for frame in range(no_frames):
                print("Time: {}   Average: {}   Standard deviation: {}".format(time_list[frame], np.rad2deg(avg_molecule_angles[frame]), std_devs[frame]), file=text_file)
    else:
        np.save(args.output, [time_list, np.rad2deg(avg_molecule_angles), std_devs])
    
