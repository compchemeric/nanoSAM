import mdtraj as md
import numpy as np
from scipy.spatial import distance
import argparse


def lattice_constants(trajectory):
    thio_indices = trajectory.topology.select('name THIO or name S')
    
    averaged_six_nearest_neighbours= []
    std_dev = []
    
    for frame in range(trajectory.n_frames):
        thio_x = trajectory.xyz[frame,thio_indices,0]
        thio_y = trajectory.xyz[frame,thio_indices,1]
        box_length_x = trajectory.unitcell_lengths[frame,0]
        box_length_y = trajectory.unitcell_lengths[frame,1]
        
        #Reshape to 2D-Array first because pdist requires it
        dx=distance.pdist(thio_x[:,None])
        dy=distance.pdist(thio_y[:,None])
        
        #Check for wrap-around nearest neighbours
        dx[dx>.5*box_length_x]=abs(dx[dx>.5*box_length_x]-box_length_x)
        dy[dy>=.5*box_length_y]=abs(dy[dy>=.5*box_length_y]-box_length_y) 
            
        d=distance.squareform(np.sqrt(dx**2+dy**2))
        
        #Using partition to get the lowest 7 distances for each thio, then sort to get everything but the smallest distance
        six_nearest_neighbours = np.sort(np.partition(d, 7, axis=1)[:,0:7], axis=1)[:,1:7]
        
        averaged_six_nearest_neighbours.append(np.nanmean(six_nearest_neighbours))
        std_dev.append(np.std(six_nearest_neighbours))
                                               
    return averaged_six_nearest_neighbours, std_dev
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Lattice constant calculation')
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
    
    args = parser.parse_args()
    
    averaged_six_nearest_neighbours = []
    std_dev = []
    no_frames = 0
    time = []
    
    #Handle streaming or completely loading the file
    if args.stream:
        for chunk in md.iterload(args.trajectory,top=args.topology, chunk = args.chunks, stride=args.frameskip):
            averaged_six_nearest_neighbours_chunk, std_dev_chunk = lattice_constants(chunk)
            averaged_six_nearest_neighbours.extend(averaged_six_nearest_neighbours_chunk)
            std_dev.extend(std_dev_chunk)
            time.extend(chunk.time)
            no_frames += chunk.n_frames
        
    else:
        trajectory = md.load(args.trajectory,top=args.topology, stride=args.frameskip)
        
        time = trajectory.time
        no_frames = trajectory.n_frames
        
        averaged_six_nearest_neighbours, std_dev = lattice_constants(trajectory)
    
    #Print the average
    print(np.mean(averaged_six_nearest_neighbours))
    
    #Save the results
    if args.format != 'numpy':
        with open(args.output, "w") as text_file:
            for frame in range(no_frames):
                print("Time: {}   Average: {}   Standard deviation: {}".format(time[frame], averaged_six_nearest_neighbours[frame], std_dev[frame]), file=text_file)
    else:
        np.save(args.output, [time, averaged_six_nearest_neighbours, std_dev])