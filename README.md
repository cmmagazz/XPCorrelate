 XPCorrelate
#
Code to register and analyse multidimensional maps e.g. nanoindentation maps vs EBSD, etc.
The default setup is to analyse Nanoindentation mapping against EBSD data, and if desired, EPMA data.
Data formats are expected as:
    - Nanoindentation: .mat file following the format given in XPImport. See below. 
    - EBSD: ctf or h5 file. 
    - EPMA: currently using a .tiff image for this due to unit difficulties.

Requirements, external: 
- Matlab
- MTEX v5.0.3 
- For second correction, ECC @ https://github.com/harshagurnani/ECC-for-Image-registration
- For CTF saving in the second correction, credit to Dr. Azdiar A. Gazder and Dr. Frank Niessen for https://github.com/frankNiessen/exportCTF

Requirements, internal: 
- XPImport @ https://github.com/cmmagazz/XPressImport if you haven't already prepped the nanoindentation data

#
How to use XPCorrelate
Most of the process flow is described in XPC_InputDeck.m 
This starts with several sections: 
1) Clearing workspace and loading relevant packages (MTEX and ECC etc)
2) User input section. Largely self explanatory in place and short enough to read in place
3) Loading and cleaning: loads the nanoindentation mapping and ebsd data. Formats the filepath for saving things properly, and crops nanoindentation data if need be. 
4) Run the EBSD registration script. See this for details. Uses affine transformation based on User input points (>=4)
   This creates the variable *datastack*, a struct where each variable should be an array of size no of x indents * no of y indents. 
   The convention used is that data should be structured in grids which start at 0,0, and X(2,1) is larger than X(1,1), and Y(1,2) is larger than Y(1,1).
   Datastack structure: 
 datastack.X               %X position
 datastack.Y               %Y position
 datastack.S               %Surface displacement
 datastack.D               %Depth
 datastack.L               %Load
 datastack.M               %Modulus
 datastack.St              %Stiffness^2/Load
 datastack.H               %Hardness
 datastack.phi1            %phi1  (ebsd)
 datastack.Phi             %Phi   (ebsd)
 datastack.phi2            %phi2  (ebsd)
 datastack.phase           %phase (ebsd)
 datastack.BCebsd          %Band Contrast (ebsd)
5) Save figures
6) Run step 4 and 5 for EPMA data. 
7) Save the registered ebsd data
8) Grain boundary analysis: creates some new 
 datastack.GBD             %Distance to nearest grain boundary
 datastack.GBSZ            %Size of grain in which point resides
 datastack.gID             %ID of grain in which point resides
9) Save all variables to .mat file. 

#
Other scripts
Naming convention: 
- f_: function. Primarily related to the main XPC_InputDeck, usually run as is. 
- g_: Second correction. Comes after f. 
- l_: Least square function, to be edited and used post-manipulation.
- p_: post-correction data manipulators. 
- v_: visualisation. 
- z_: Scrap, but handy. 

From here on, scripts are run ad-hoc as needed. In order as needed to recreate the results in CM Magazzeni et al 2020:

PCDM I.e. Post-Corrected Data Manipulator
    Needed to plot various graphs, and contains other snippets which weren't used in publication but are fun or useful. 

lsqfitting-XXX
    lsqfitting.m : the first least square fitting script. Used to plot H vs Phi, with or without filtering based on GBD. 
    lsqfittingOXYGEN_XXX.m : various least square fitting scripts for O vs H. 
