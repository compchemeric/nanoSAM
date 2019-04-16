import mdtraj as md
import numpy as np
import argparse
from scipy.spatial import distance
from scipy.cluster.hierarchy import fcluster, linkage


def closest_node(node, nodes, x_len, y_len):
    #Emulate cyclic repetition
    nodes_x = nodes[:,0]
    nodes_x[nodes_x-node[0] > 0.5*x_len] -= 0.5*x_len
    nodes_x[node[0]-nodes_x > 0.5*x_len] +=  0.5*x_len
    
    nodes_y = nodes[:,1]
    nodes_y[nodes_y-node[1] > 0.5*y_len] -= 0.5*y_len
    nodes_y[node[1]-nodes_y > 0.5*y_len] += 0.5*y_len
    
    merged_nodes = np.array((nodes_x,nodes_y)).T
    
    #Get index of min distance
    return distance.cdist([node], merged_nodes).argmin()


def voronoi(coordinates, x_len, y_len, bins):
    x_dimension = np.linspace(0, x_len, bins);
    y_dimension = np.linspace(0, y_len, bins);
    
    vertex = np.empty(shape=(bins,bins), dtype=np.int8)
    
    for i_x in range(bins):
        for i_y in range(bins):
            vertex[i_x][i_y] = closest_node([x_dimension[i_x],y_dimension[i_y]], coordinates, x_len, y_len)
            
    return vertex


def cluster(atoms):
    distances = distance.pdist(atoms - np.mean(atoms,axis=0))
    clusters = fcluster(linkage(distances,  method='single'), 2, criterion='maxclust')
 
    return atoms[clusters==1], atoms[clusters==2]
    

def get_thickness_map(grid, cluster, other_cluster, x_len, y_len):
    thickness_grid = np.empty(shape=grid.shape)
    
    for i_x in range(len(grid)):
        for i_y in range(len(grid)):
            nearest_po_in_other_layer = closest_node(cluster[grid[i_x][i_y], 0:2], other_cluster[:, 0:2], x_len, y_len)
            thickness_grid[i_x][i_y] = abs(cluster[grid[i_x][i_y], 2] - other_cluster[nearest_po_in_other_layer, 2])
    
    return thickness_grid
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Bilayer thickness')
    parser.add_argument('--trajectory', help='Trajectory input path', required=True)
    parser.add_argument('--topology', help='Topology input path', required=True)
    parser.add_argument('--output', help='Path to output', required=True)
    
    parser.add_argument('--txt', dest='format', help='Output as txt', action='store_const', const='txt')
    parser.add_argument('--numpy', dest='format', help='Output as numpy readable file', action='store_const', const='numpy')
    parser.set_defaults(format='txt')
    
    parser.add_argument('--frameskip', dest='frameskip', help='Only evaluate every nth frame', type=int)
    parser.set_defaults(frameskip=1)
    
    parser.add_argument('--bins', dest='bins', help='Number of bins for voronoi', type=int)
    parser.set_defaults(bins = 100)
    
    args = parser.parse_args()
    
    trajectory = md.load(args.trajectory,top=args.topology, stride=args.frameskip)
    
    thickness_map_list = []
    avg_thickness_list = []
    
    for frame_num in range(trajectory.n_frames):
        
        box_length_x = trajectory.unitcell_lengths[frame_num,0]
        box_length_y = trajectory.unitcell_lengths[frame_num,1]
    
        po_xyz = trajectory.xyz[frame_num, trajectory.topology.select('name PO4')]
    
        cluster1, cluster2 = cluster(po_xyz)
        
        #Use both layers for thickness calculations and average it
        grid1 = voronoi(cluster1[:, 0:2], box_length_x, box_length_y, args.bins)
        thickness1 = get_thickness_map(grid1, cluster1, cluster2, box_length_x, box_length_y)
        
        grid2 = voronoi(cluster2[:, 0:2], box_length_x, box_length_y, args.bins)
        thickness2 = get_thickness_map(grid2, cluster2, cluster1, box_length_x, box_length_y)
    
        average_thickness_map = np.mean([thickness1, thickness2], axis=0)
        
        thickness_map_list.append(average_thickness_map)
        avg_thickness_list.append(np.mean(average_thickness_map))
    
    average_thickness_map_over_all_frames = np.mean(thickness_map_list, axis=0)
    
    #Save the results
    if args.format != 'numpy':
        with open(args.output, "w") as text_file:
            for frame in range(trajectory.n_frames):
                print("Time: {}   Average: {}".format(trajectory.time[frame], avg_thickness_list[frame]), file=text_file)
                
            print("----------------------------------------------------------------------------\n", file=text_file)
            
            for i_x in range(len(average_thickness_map_over_all_frames)):
                for i_y in range(len(average_thickness_map_over_all_frames)):
                    print("X: {}   Y: {}   Thickness: {}".format(i_x, i_y, average_thickness_map_over_all_frames[i_x][i_y]), file=text_file)
                    
    else:
        np.save(args.output, [trajectory.time, avg_thickness_list, average_thickness_map_over_all_frames])