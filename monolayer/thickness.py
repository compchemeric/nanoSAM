import numpy as np
import argparse


#Creates a line from two points
def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C


#Returns x and y coordinates of a intersection of 2 lines
def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False
    
    
def thickness(hist1, hist2, z_length):
    hist = hist1 - hist2
    bins = len(hist)
    
    sign_changes = []
    intersections = []
    
    #Determine all possible locations for intersections
    for i in range(bins-2):
        if np.sign(hist[i]) != np.sign(hist[i+1]):
            sign_changes.append(i)
    
    #Linear interpolation for more accurate intersection point
    for intersect in sign_changes:
        line1 = line((intersect, hist1[intersect]),(intersect+1, hist1[intersect+1]))
        line2 = line((intersect, hist2[intersect]),(intersect+1, hist2[intersect+1]))
    
        intersections.append(intersection(line1, line2))
    
    #Calculating possible thicknesses
    possible_thickness = []
    
    for i in range(len(intersections)):
        for j in range(i+1, len(intersections)):
            thickness = (intersections[j][0]-intersections[i][0]) * z_length/bins
    
            possible_thickness.append(thickness)
    
    return possible_thickness


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Thickness calculation')
    parser.add_argument('--input', help='input path', required=True)
    
    parser.add_argument('--histograms', help='List of historgrams to be considered, like so: (a,b)(c,d)')
    
    args = parser.parse_args()
    
    histograms = np.loadtxt(args.input)
    
    #Parsing
    z_length = histograms[-1][0]
    histograms = histograms[:-1]
     
    histogram_list = args.histograms
    histogram_list = histogram_list.replace(" ", "").replace(")(", "#").replace(")", "").replace("(", "").split("#")
    
    #Looping over all Tuples
    counter = 1
    for tuple_of_histograms in histogram_list:
        tuple_of_histograms = tuple_of_histograms.split(",")
        
        possible_thickness = thickness(histograms[int(tuple_of_histograms[0])], histograms[int(tuple_of_histograms[1])], z_length)
        
        print("Possible thickness for tuple {}:".format(counter))
        
        for t in possible_thickness:
            print(t)
        
   
