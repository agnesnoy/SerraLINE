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

SerraLINE is a software for the calculation of bending angles at 
different lengths (in bp steps)) from the global molecular contour 
of simulations of DNA. This molecular contour filters out local 
irregularities caused by the helical periodicity of DNA. 
For instance, the molecular contours were calculated by WrLINE, 
but any contour with the required format could be processed 
by SerraLINE. 

Closed (circular) and opened (linear) DNAs can be processed by 
SerraLINE, but the user is required to specify this in the
input file SerraLINE.in. 

Bending angles are obtained through tangent vectors, where the 
tangent vector of the last base-pair is lost for opened structures.
This causes to have (N-1)(N-2)/2 bendings for opened structures, 
whereas N(N-1) for closed structures, where the tangent vectors of
all base-pairs can be obtained. 

SerraLINE is composed of two main programs:

  - SerraLINE: Main process that calculates bending angles at 
               different lengths (in bp steps) and outputs 
               averages and standard deviations of these angles.

  - Extract: Supportive process which filters out angles at a 
             particular length of interest (indicated by the user) 
             and sorts out the positions according to the 
             hinge point.
           
Additionally to the standard method, SerraLINE can project the 
input structure to the plane that best fits the molecule or a 
given region specified by the user.
This projection method mimics experiments where structures are 
visualized on a two dimensional plane, and quantities such as width, 
height and aspect ratio (width/height) can be calculated.

For both methods, bending angles are calculated with the
same criteria.

An example is provided for running SerraLINE and processing 
results.

___________________________________________________________________




__________________________________________________________________
                       INPUTS
                       ------

For running SerraLINE you only need a trajectory file:


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


For running the main process, the user has to indicate the path
to the data files as well as the type of structure in the next
input file:

  - SerraLINE.in .
    Seven inputs need to be indicated in this file. 
    1.- The user first needs to specify the type of structure, 
        typing 0 for opened structure (linear DNA) and 1 for
        closed structure (circular DNA).
    2.- In case that you have a topology file corresponding to 
        a single stranded DNA, then type 1, whereas if the
        topology corresponds to a double stranded DNA, then type
        2. If you don't have a topology file, then type 0.
    3.- In case that you type 0 in the previous input, then 
        you are required to indicate the number of base-pairs in
        the trajectory file. With this, SerraLINE can determine
        the number of frames in the simulation. If you didn't
        type 0 in the previous input, then SerraLINE will ignore
        this input.
    4.- Method selection. Type 1 if you want to perform the
        standard method. For the projection method, by typing
        0 the whole structure will be projected to the plane 
        that best fits, or, specifying with ' ', you can 
        fit a plane to a selection of base-pairs. This 
        selection of points should be at least 3 bps.
    5.- Path to the topology file. This input will be ignored if
        you indicated that you don't have a topology.
    6.- Path to the trajectory file. This trajectory can be in
        amber format (*crd or *x) and WrLINE outputs format
        ( *3col and *xyz).
    7.- In case the projection method was selected, the projected
        trajectory can be written in xyz or crd format. Type 1
        for writing the projected trajectory or 0 for not.

Another input file is required for executing the Extract process:

  - extract.in .
    Only two inputs are required this time.
    1.- Path to SerraLINE's output file (SerraLINE.out)
    2.- Length of interest (in bp steps) to filter angles.


These input files (*in) can be renamed.

___________________________________________________________________




___________________________________________________________________
                       EXECUTION
                       ---------

To execute SerraLINE, simply type: SerraLINE < SerraLINE.in

This will generate the output "SerraLINE.out"

If you want to extract lengths (for processing data), then 
specify the length of interest in extract.in as well as the path to
SerraLINE's output file (SerraLINE.out)
Then type: ./Extract < extract.in

This will filter out parameters at yout length of interest and will
also sort the data (it will be ready to plot!).

___________________________________________________________________




___________________________________________________________________
                       OUTPUTS
                       -------

The main process produces the next output:

  - SerraLINE.out.
    Contains information of the structure analysed and average and 
    standard deviations of bending angles at different lengths.
    If the projection method was used, average and standard 
    deviation of width, height and aspect ratio will be given at
    the top.

And the extraction process (Extract) produces the following output:

  - subfragment_$l.out .
    Bending angles at length of interest $l (indicated in bp steps)
    and rearrenged according to the hinge positions.

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

You can eather modify extract.in and extract bending angles at a
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



