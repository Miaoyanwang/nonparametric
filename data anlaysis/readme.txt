###########################
bbnet68_spatial_orientation.RData
Reference: Dunson_AOAS19.pdf.
Compare your analysis result to two papers, named Dunson_*19.pdf, in the reading folder.

A: a 68 × 68 × 212 binary tensor consisting of structural connectivity patterns among 68 brain regions for 212 individuals from Human Connectome Project (HCP). The tensor entries encode the presence or absence of fiber connections between those 68 brain regions for each of the 212 individuals.

VSPLOT (Variable Short Penn Line Orientation Test): two line segments are presented on the screen and participants are asked to rotate a movable line so that it is parallel to the fixed line. The sample consists of 106 individuals with high (top 10%) VSPLOT scores and 106 with low (bottom 10%) VSPLOT scores.

Node file for VSPLOT: ROI_68_lobe_node_renamed.txt

Node file for IQ data: node file ``Desikan-Killiany68.node''


###########################
BrainNet is a useful software for visualizing brain connections.

Instructions (The detailed instruction is brainNet_Manual.pdf):
1. Open matlab and type "BrainNet"
2. Load files:
Surface file: BrainMesh_ICBM152.nv

Data file (node): depending on the input file. 

Data file (edge): A user-specified input. This is a matrix file representing the brain connection. The edge can be raw connection strength, estimated strength, truncated strength, or any measures you want to plot over the brain edges.

Note: (You may need enable "all files *" option by clicking the "options" in the pop-up console.)
