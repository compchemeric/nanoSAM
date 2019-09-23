import mdtraj as md
import numpy as np
import argparse


def set_amu_cg():
    return { "W": 24, "WM": 24, "WP": 24, "THIO": 72, "C1A": 72, "C2A": 72, "C3A": 72, "C4A": 72, "ACA": 72, 
            "ACB": 72, "C1B": 72, "C2B": 72, "C3B": 72, "PE1": 45, "PE2": 45, "PE3": 45, "PE4": 45, "PE5": 45, "PE6": 45,
            #For residues, this is the average amu of the residue, so the sum of the amus of all atoms in the residue divided by the number of atoms in the residue
            "TE8": 846/14, "TE16": 990/16, "FIL": 405/16, "PW": 72/3,
    }
    

def set_amu_aa():
    return { 'H': 1.0, 'C':12.0, 'O': 16.0, 'N': 14.0, 'S':32.0, 'OW': 16.0, 'HW': 1.0,
            #For residues, this is the average amu of the residue, so the sum of the amus of all atoms in the residue divided by the number of atoms in the residue
            'FIL': 331/59, 'TE8': 592/117, 'TE16':0, 'HOH': 18/3
    }
    

def histogram_of_atoms(trajectory, bins, atoms_to_select, residues_to_select, pdf, amu):
    selections = []
    
    #Parse atom selection from command line
    if atoms_to_select:
        atoms_to_select = atoms_to_select.upper()
        atoms_to_select = atoms_to_select.replace(" ", "")
        atoms_to_select = atoms_to_select.split(',')
    
        #Execute selection
        for atom in atoms_to_select:
            select_string = 'name ' + atom
            sel = trajectory.topology.select(select_string)
            if sel.any():
                selections.append((sel, atom))
    
    #Parse residue selection from command line
    if residues_to_select:
        residues_to_select = residues_to_select.upper()
        residues_to_select = residues_to_select.replace(" ", "")
        residues_to_select = residues_to_select.split(',')
        
        #Execute selection
        for residue in residues_to_select:
            select_string = 'resname ' + residue
            sel = trajectory.topology.select(select_string)
            if sel.any():
                selections.append((sel, residue))
    
    #Exit if no selections are made
    if not selections:
        raise Exception('No selections can be found')
        return
    
    max_box_length_x = np.max(trajectory.unitcell_lengths[:,0])
    max_box_length_y = np.max(trajectory.unitcell_lengths[:,1])
    max_box_length_z = np.max(trajectory.unitcell_lengths[:,2])
    
    volume_of_bin_in_m3 = max_box_length_x*max_box_length_y*max_box_length_z/bins * 10 ** -27
    amu_to_kg = 1.66054 * 10 ** -27
    
    histograms = []
        
    for selection in selections:
        z = trajectory.xyz[:,selection[0],2]
        
        #Create a pdf
        if pdf:
            hist, _ = np.histogram(z, bins, density=True, range=(0, max_box_length_z))
            hist = hist.astype(float)/sum(hist)
            histograms.append(hist)
        #Create a Kg/m³ histogram
        else:
            hist, _ = np.histogram(z, bins, density=False, range=(0, max_box_length_z))
            
            #Multiply by the amu given in the dictionary or raise error if amu can not be found
            element = ''.join(filter(lambda c: not c.isdigit(), selection[1]))
            if selection[1] in amu:
                hist = hist * amu[selection[1]] 
            elif element in amu:
                hist = hist * amu[element]
            else:
                raise Exception('Amu for atom can not be found')

            hist = hist * amu_to_kg / volume_of_bin_in_m3 #Convert to Kg/m³
            hist = hist / trajectory.n_frames
            histograms.append(hist)

    return histograms


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Histogram calculation')
    parser.add_argument('--trajectory', help='Trajectory input path', required=True)
    parser.add_argument('--topology', help='Topology input path', required=True)
    parser.add_argument('--output', help='Path to output', required=True)
    
    parser.add_argument('--frameskip', dest='frameskip', help='Only evaluate every nth frame', type=int)
    parser.set_defaults(frameskip=1)
    
    parser.add_argument('--bins', dest='bins', help='Number of bins to use', type=int)
    parser.set_defaults(bins=10)
    
    parser.add_argument('--atoms', help='List of atoms to be considered, comma separated')
    parser.add_argument('--residues', help='List of residues to be considered, comma separated')
    
    parser.add_argument('--pdf', dest='pdf', help='Create a probability density function instead of a histogram with atom mass', action='store_true')
    parser.set_defaults(pdf=False)
    
     parser.add_argument('--verbose', dest='verbose', help='Verbose', action='store_const', const=True)
    parser.set_defaults(verbose=False)
    
    parser.add_argument('--resolution', dest='resolution', help='Coarse grained(CG) or AA')
    parser.set_defaults(resolution='cg')
    
    args = parser.parse_args()
    
    histograms = []
    
    if args.verbose:
        print("Loading trajectory")
        
    trajectory = md.load(args.trajectory,top=args.topology, stride=args.frameskip)
    
    if args.verbose:
        print("Trajectory loaded")
        
    #Parse resolution
    amu = {}
    
    if args.resolution.upper() == "CG":
        amu = set_amu_cg()
    elif args.resolution.upper() == "AA":
        amu = set_amu_aa()
    else:
        raise Exception('Resolution not valid')
    
    #Execute
    histograms = histogram_of_atoms(trajectory, args.bins, args.atoms, args.residues, args.pdf, amu)
    
    #Save the results
    last_line = [np.zeros(len(histograms[0]))]
    last_line[0][0] = np.max(trajectory.unitcell_lengths[:,2])
    
    with open(args.output, "w") as file:
        np.savetxt(file, histograms)
        np.savetxt(file, last_line)
    
    if args.verbose:
        print("File saved to: {}", args.output)