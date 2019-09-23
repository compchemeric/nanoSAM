# nanoSAM
Nanoscale Characterization of Self-Assembled-Monolayers and tethered Lipid Bilayers and Vesicles

The scripts above were used in our article "Simulation of Mixed Self-Assembled Monolayers on Gold: Effect of Terminal Alkyl Anchor Chain and Monolayer Composition", published in J. Phys. Chem. B, 2018, 122 (31), pp 7699–7710 under the DOI: 0.1021/acs.jpcb.8b05075. 

They should work with any molecular dynamics trajcetory data of self-assembled monolayer systems and allow for a quick computation of lattice constant, tilt angles, and thickness, order paramters as well as trans- and helical ratios.

The output can be written to either text or numpy files. Examples and some more documentation will follow soon.  

##General

Short explanation of the command line parameters:

    --topology Topology input file

    --trajectory Trajectory input file

    --output File the output will be written to

    --stream (When available) Enable streaming of the trajectory input file, instead of loading it all at once, reducing the ram-usage.

    --chunks Size of the chunks to be loaded when streaming is enabled, in frames.

    --frameskip Only load every n-th frame of the trajectory.

    --verbose Display additional information helpful for debugging.

    --txt Output as humanly readable text.

    --numpy Output as native numpy file for easy loading for further processing.

##Monolayers

##Bilayers
####thickness.py
Create a map of a bilayers thickness (P to P atom) by creating a 2-D grid and assigning every gridpoint the nearest calculated thickness.

    --bins: Number of points in both directions of the grid. A higher number of points means a higher accuracy, but the more points, the slower the computation.

Example Output as textfile: 
X: 1   Y: 1  Thickness: 2
X: 1   Y: 2  Thickness: 3
X: 2   Y: 1  Thickness: 1
X: 2   Y: 2  Thickness: 3

Output as a numpy file is recommended for easy visualtization with matplotlib.