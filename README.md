___________________________________________________________________
                        \       _
			'       '
                \       ---------        /
		|--------SERRALINE-------|
                        ---------
                        |       |
                        -       -
___________________________________________________________________
___________________________________________________________________




___________________________________________________________________
                       REQUIREMENTS
                       ------------

You only need a fortran compiler (gfortran), that's it.

___________________________________________________________________




___________________________________________________________________
                       COMPILATION
                       -----------

Type make SerraLINE for compiling SerraLINE.
     make Extract for compiling the extraction process.
     or make all for compiling both.
___________________________________________________________________




___________________________________________________________________
                       OVERVIEW
                       --------

SerraLINE is a software that calculates bending angles, width, height,
aspect ratio and deviatin from planarity of DNA molecules using the 
global molecular contour WrLINE. 

The molecular contour WrLINE defines a point for each bp.

SerraLINE can process closed (circular) and opened (linear) DNAs.

Bending angles are measured between two tangent vectors that can be separated 
by different number of bps (or points).

Tangent vectors are constructed by the vectors joining two points that 
can be consecutive or arbitrarily separated.
This latter option is included to further approach microscopy imaging and 
analysis with a measurable limited resolution. 

SerraLINE can project the WrLINE contour to the plane that best fits the 
molecule or to the given region specified by the user.
This projection method mimics experiments where structures are 
visualized on a two dimensional plane.
Then, quantities such as width, height, aspect ratio (width/height)
and deviation from planarity can be calculated.

For both options (projection or not), bending angles are calculated 
following the same criteria.

SerraLINE is composed of two main programs:

  - SerraLINE: Main program that calculates bending angles at 
               different lengths (in bp steps) and outputs 
               averages and standard deviations of these angles.
               When the projection is actived, then it outputs 
               width, height, aspect ratio and deviation from 
               planarity
               

  - Extract: Supportive process which filters out angles at a 
             particular length of interest (indicated by the user) 
             and sorts out the positions according to the 
             middle point.
           
An example is provided for running SerraLINE and processing 
results.

___________________________________________________________________




__________________________________________________________________
                       SerraLINE
                       --------

For running SerraLINE you need:


  - Trajectory file (Essential): 
    A trajectory file of the molecular contour. 
    Three formats are supported by SerraLINE. Either *.xyz and 
    *.3col which basically are outputs of WrLINE and *.crd (or *.x) 
    which is the trajectory in amber style format (10F8.3). 


  - Topology file (Optional).
    A topology file in amber style format. By using a topology
    file SerraLINE reads the DNA sequence and determines the number 
    of bases from this file.
    SerraLINE can also work without the topology file, but the 
    sequence will not be obtained and the number of base-pairs will
    need to be provided by the user (in SerraLINE.in).

  - SerraLINE.in 
    Seven keys need to be indicated in this file. 
    The different options are also specificed there.
    1.- The user first needs to specify the type of structure, 
        typing 0 for opened structure (linear DNA) and 1 for
        closed structure (circular DNA).
    2.- In case you have a topology file corresponding to 
        a single stranded DNA, type 1.
        If the topology corresponds to a double stranded DNA, then type
        2. 
        If you don't have a topology file, then type 0.
    3.- In case that you type 0 in the previous input, then 
        you are required to indicate the number of base-pairs
        of your molecule. With this, SerraLINE can determine
        the number of frames in the simulation. 
        If you didn't type 0 in the previous input, 
        then SerraLINE will ignore this input.
    4.- Projection of the trajectory on to the best fitted plane. 
        Type 1 if you do not want to project.
        Type 0 if you want to project
        Type 'x:y' if you want to fit the plane to a selection o bp *x and y) 
        This selection of points should be at least 3 bps.
    5.- Tangent length specification. Type the length l that
        defines how the tangent vectors are constructed.
    6.- Path to the topology file. This input will be ignored if
        you indicated that you don't have a topology.
    7.- Path to the trajectory file. This trajectory can be in
        amber format (*crd or *x) and WrLINE outputs format
        ( *3col and *xyz).
    8.- In case the projection method was selected, the projected
        trajectory can be written in xyz or crd format. Type 1
        for writing the projected trajectory or 0 for not.

To execute SerraLINE, simply type: 
./SerraLINE < SerraLINE.in

The SerraLINE produces the next output:

  - SerraLINE.out.
    Contains information of bending angles at different lengths.
    If the projection method is used, width, height, aspect ratio 
    and deviation from planarity are given at the top.
    It calculates such quanitities for each frame and outputs average 
    standard deviation.
    For deviation from planarity, the program prints the distance 
    from the plane averaged along the molecule and the point that 
    is farther apart. It outputs the absolute distance (in angstroms)
    and the relative distance (in %) with height

___________________________________________________________________



___________________________________________________________________
                       Extract
                       ---------

This utility program extracts the bend profile along the polymer 
at a particular length. It has the following input

  - extract.in .
    Only two inputs are required this time.
    1.- Path to SerraLINE's output file (SerraLINE.out)
    2.- Length of the bendinterest (in bp steps) of the ben.


Then type: 
./Extract < extract.in

This will filter out bend angles at yout length of interest and make it ready to plot 
on the following output:

  - subfragment_$l.out .
    Bending angles at length of interest $l (indicated in bp steps)
    and rearrenged according to the middle positions.

___________________________________________________________________




___________________________________________________________________
                       EXAMPLE
                       -------

To help the user familiarise with SerraLINE, we use outputs from 
WrLINE's example:

 -C1.xyz and C1.3col contain the trajectory, you can use any of 
  those specifying which one in SerraLINE.in

 -test.prmtop is the topology file.

In SerraLINE.in, it already specifies the type of structure (circular),
no topology is being provided, the number of base-pairs and the
path to the trajectory file. This file performs the standard method 
(without projection).

The next process will compile everything and run the main process

Compilation: make all 
Execution: ./SerraLINE < SerraLINE.in

This will produce SerraLINE.out, which contains bending angles
and information of the structure analysed (no sequence is provided
since in SerraLINE.in it is indicated that we won't use the topology
file).

You can either modify extract.in and extract bending angles at a
particular length or execute run-example.sh, which will run Extract
and filter angles of lengths from 1 to 5 bp steps.

We also provided a simple python script to plot the bending angles
at lengths from 1 to 5 bp steps.

You can run this example by typing:

./run-examples

and then

ipython plot_angles.py

This will produce a graph where each colour corresponds to angles
at particular lengths, and shows how strong bends arise from lower
lengths to higher lengths. The position of the strong bends, 
correspond to the U-turns in the DNA.

___________________________________________________________________



