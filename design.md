# Introduction

In this directory, /lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS, I want to build a script routine to sample rotamer states of protein residues lining a binding site cavity using Python, pyKVfinder package, and the Schrodinger Maestro tool suite. The goal is to take a reference ligand-bound structure from the PDB, which will be provided by the user, along with a series of all-atom structural samples of the same protein generated using the BioEmu model. For all generated samples, we want to identify the amino acid residues lining the binding site in the sample by first using the reference ligand-bound structure to select the bound ligand, select residues closer than 7 angstrom distance from the ligand, and use that list of residues to define the binding site cavity in the sample structures. (The sample structures do not have a ligand bound, so we need to use a reference experimental ligand-bound structure.) Once we have a list of the residues lining the binding cavity from the experimental PDB structure, we can focus on modifying the sample structures.

I use "amino acid" and "residue" interchangeably in this document. "Binding site", "cavity", "pocket", "ligand site", and such terms are also used interchangeably. "Context-dependent rotamer library" and "backbone-dependent rotamer library" both refer to the same thing, which is the schrodinger.protein.rotamers.Rotamers object/module. This module takes two arguments, the structure (in Maestro format, from schrodinger.structure.StructureReader.read(file) ) and one atom object from the residue you want to load the corresponding rotamer library for.

The goal of the scripts will be to maximize the measured volume of the cavity in each sample structure that contacts all of the identified binding-site-lining residues from the experimental structure. This maximization algorithm should proceed in two phases. We will do this by using Schrodinger's Maestro suite to first mutate all of the residues in the binding site out to Alanine, one by one, and measure the effect that this mutation has on the cavity volume. In doing this, we want to create a sorted list of the binding site residues - sorted by their impact on the cavity volume when mutated. The single-residue mutations that lead to the LARGEST positive changes in cavity volume (i.e., those that maximally INCREASE the volume of the cavity, as identified by pyKVfinder) should be the ones prioritized for rotamer sampling in the second phase of the algorithm. The single-residue mutations that lead to NO changes in cavity volume (i.e., those that do not impact the pyKVfinder measurement of cavity volume upon mutation to Ala) should be deleted from the residue list, as deleting the sidechain from that residue leads to no steric benefit for the binding pocket.

The first phase of the algorithm is the mutation and residue sorting step specified above. The second phase of the algorithm will involve rotamer sampling again using the Maestro suite. Proceeding down the list of sorted residues one by one (i.e., starting with the residue whose mutation to Alanine maximally increases the pocket volume with respect to all other residue mutations), we want to go through Maestro's pre-built context-dependent rotamer state library for that particular residue, apply all the states, and create a list of rotamer states for EACH residue in the binding-site list; this list of per-residue rotamer states should again be sorted, in decreasing order, by the change in cavity volume that each rotamer state enables. I.e., the first rotamer state in this sorted list (for each residue) should be the rotamer state of the residue that maximally increases the cavity volume (from a single sidechain rotamer transition), and so on. While we do this, we also want to extract the frequency of observing each rotamer state in the above list for that particular residue from Maestro's context-dependent rotamer library. Sort this list of rotamer state frequencies according to the same criteria as above, i.e., sort the rotamer state frequency list by the effect that that rotamer state has on cavity volume, in decreasing order. In this way, each item at each index of the list of rotamer states, as well as the list of rotamer state frequencies, should correspond to the same rotamer state.

At the end of the algorithm, we should take the rotamer states for each cavity-lining residue that maximize pocket volume, and we should apply all of those states to the sample structure to maximize the cavity volume. When this is done, please run pyKVfinder one last time on the modified sample structure to ENSURE that the cavity volume is LARGER than the original sample structure; please report in cubic angstroms how much the cavity volume has been enlarged. In addition, we should make sure that the final modified sample structure only contains ONE identified pocket (identified by pyKVfinder) that is surrounded by the pocket-lining residues. Sometimes, one large binding pocket/cavity is recognized as two smaller distinct cavities because of a sidechain state problem that sterically splits the cavity up. When this occurs, which you will see in the example sample structure (s0067_rec.pdb), we want to identify ALL cavities, no matter how small, that are in contact with at least two of the residues identified to line the binding pocket from the experimental reference structure. The effective "cavity volume" to maximize is the SUM of the volumes of all of these residue-contacting cavities. We know we have reached a good modified structural sample when only ONE residue-contacting cavity exists, and the volume of that cavity is larger than the initial cavity volume (whether it is a sum of distinct cavities in the original sample structure or otherwise).

# Base considerations for software environment

We will likely create python scripts following the structure and syntax in the files found in /lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/example_script. Before testing or running any of these Python files, however, we must be sure to activate the Schrodinger license environment by running these two commands:

export SB_BASE_OVERRIDE=/lustre/fs4/ruit/scratch/rebecca/LYU_SBGRID/INST
source /lustre/fs4/ruit/scratch/rebecca/LYU_SBGRID/INST/sbgrid.shrc

In addition, please always activate and use the conda environment called "rotamer_dms", which has pyKVFinder installed.

Once the environment has been sourced, we can run any relevant script with:

$SCHRODINGER/run python3 script.py

If you want to set up a shell script that makes running this script far easier, and also lets the user define any relevant parameters (e.g. location of experimental reference structure PDB file, location of sample structure PDB file, threshold of number of top residues to select for rotamer sampling), please do so. It would help a lot. 

# Design choices and workflow

There are three core principles in the design of this algorithm.

The first: "The optimal solution (which is an all-atom protein structure in PDB format) is the one which jointly maximizes pocket volume while minimizing the change in conformational free energy." We will not add a conformational free energy consideration in this algorithm yet; I will do that manually. However, we want to maximize pocket volume.

The second: "The residues [lining the binding site in the experimental ligand-bound structure] leading to the largest positive volume deviation upon mutation to Alanine should be the residues considered for rotamer transitions." Please give the user an option to select a percentage threshold of top-volume-affecting residues in the binding pocket to modify with rotamer state transitions in the second phase of the algorithm.

The third: "This algorithm should run as fast as possible. Calls to the pyKVfinder package for cavity detection and volume measurement may incur substantial overhead, so use this only when it is absolutely necessary to obtain a volume measurement for optimization." 

# Precise specifications for algorithm development

When we compute the volume of the cavity in any sample structure, we want to be sure that the cavities we identify with pyKVfinder correspond to the cavity in which the ligand actually binds, as observed in the experimental reference structure. pyKVfinder will likely find many cavities in any of our sample structures, and most of them will not correspond to the appropriate ligand-binding site. There is no correlation between pocket volume, depth, or any other measurement with the "correctness" of that pocket for ligand binding. We can only determine the cavity to be the "correct" one by ensuring that the residues in contact with that cavity are the same residues that are in contact with the cavity in the experimental reference structure. You will find that sometimes pyKVfinder finds multiple (two or more) cavities that contact some subset of the residues found to border the binding cavity in the experimental reference structure. If this is the case, make sure that each of those identified cavities in the sample structure contacts at least two of the residues that we know to contact the cavity in the experimental reference structure. We can then treat those separate identified cavities together; the volume to maximize is the SUM of their cavity volumes. Ideally, at the end of the algorithm, we only have ONE contiguous cavity with a volume enlarged compared to that of the initial sample structure.

Perform these steps in order:

1. Preparation steps
- Read in the experimental reference structure PDB with schrodinger.structure.StructureReader.read(file).
- Read in the sample structure PDB with schrodinger.structure.StructureReader.read(file).
- Get a list of all residues in both structures. Make sure that the two structures have IDENTICAL amino acid sequences. Throw an error if the two structures do not correspond to the same protein (i.e., if they do not have 100% sequence homology or identity).
- Select the ligand in the experimental structure. Make the user provide the residue name (three-letter identification code) of the ligand to ensure we are selecting the right thing and finding the right pocket in the experimental reference structure.
- Select all residues within 7 angstroms distance from the selected ligand in the experimental structure.
- Save the list of those residues, and their sequence indices, as our working residue list to modify/mutate in the sample structures. This is our list of residues which lines the binding cavity. Print this list of residues and name it accordingly.

2. Mutation phase of algorithm
- Run pyKVfinder's cavity detection routine on the sample structure. Search for cavities in the sample structure that have at least two residues from our working residue list in contact with them.
- If no cavities are found, return an error message indicating that no suitable cavities were found in the sample structure.
- If one or more cavities are found, print the full list of residues that contact each cavity, so I can be sure that we are looking at the right cavities. The "volume" of this binding pocket in the sample structure is the sum of the volumes of each of the distinct cavities that are found to contact the binding-pocket-lining residues. Print the initial volume (or volume sum) of the relevant binding site cavity/cavities in the sample structure in cubic angstrom, so I can verify it.
- Begin mutating each of the residues in the working residue list to Alanine, one by one. IF mutating a residue to Alanine increases the volume (or volume sum) of the binding site in the sample structure, keep it in the list, and record how much it increases the volume by (in cubic angstroms) in a separate list. IF mutating a residue to Alanine does not change the volume of the binding site in the samples structure at all, remove it from the list. Please print the name (3-letter code) and sequence index of the residue when it is removed so I know what is going on.
- With the filtered working residue list from the previous step, sort the list of residues by the impact on cavity volume it has when mutated to Ala. We will use this sorted list of residues for sidechain sampling in the second phase of the algorithm. 

Overall, the mutation phase should give us a sense of which residues within 7 angstroms of the bound ligand actually sterically affect the binding pocket. It should additionally let us isolate the residues within the list for whom cavity volume is maximally increased upon mutation to Alanine, as these are the residues that are prime candidates for rotamer sampling. 

3. Rotamer sampling phase of algorithm
- This second phase of the algorithm should operate on the modified/filtered working residue list from the first (mutation) phase of the algorithm. Because this list of residues should be sorted (in decreasing order) by their contribution to pocket volume when mutated, we can work through each residue in this filtered working list one by one.
- For every residue in this filtered working list, select this residue in the sample structure and load Maestro's context-dependent rotamer library for that particular residue. Go through every rotamer state in that library. For each rotamer state in the library loaded for each residue, apply that rotamer state to the residue, use pyKVfinder to compute the new volume (or volume sum) of the binding site cavity, and also record the frequency of observing that rotamer state if that frequency information can be extracted from Schrodinger's rotamers.Rotamers object. Create a list of rotamer states and state frequencies. Sort both lists in decreasing order of rotamer states' contributions to pocket volume, i.e., the rotamer state that leads to the largest increase in pocket volume compared to the native rotamer state in the sample structure should be the first one in the rotamer state list, and the frequency of observing that rotamer state (as extracted from Schrodinger's rotamer library) should be the first entry in the state frequency list. Make sure to record the frequency of observing the native rotamer state (that observed for the residue in consideration in the unmodified sample structure) as determined from the Schrodinger rotamer library, as well. 
- Assume the top rotamer state for each residue in consideration is the first item in the rotamer state list for that residue. Apply all such "top states" to the original sample structure. Use pyKVfinder to compute the volume of the binding site cavity as before, save the modified structure, and report the volume of the cavity in the modified structure. If pyKVfinder is still detecting more than one cavity contacting the binding site residues in the modified sample structure, save the structure anyway with a different filename - maybe append "still_multiple_pockets" to the filename - so I can inspect it.
- In its current formulation, the algorithm is done after this. Print "done", and ensure the modified structure/s are saved. If you can, while you build this workflow, it would be nice for you to save checkpoint structures and statistics (e.g. the residue lists and volume changes upon mutation or rotamer transitions) in a /checkpoints/ directory so I can validate and verify.

4. Downstream considerations: DO NOT IMPLEMENT THIS NOW.
- I want to add in functionality later to approximate the conformational free energy change from rotamer transitions (compared to the native sample structure), which is why I'm having you save rotamer state frequencies in separate lists from the Schrodinger rotamer library. I also want to add in functionality to have a user-defined potential during the rotamer sampling phase, like a Weeks-Chandler-Andersen type steric repulsion potential to ensure that chosen rotamer states do not sterically clash with each other. Don't deal with this now.

# Helpful resources

pyKVfinder Python package documentation is here:
https://lbc-lnbio.github.io/pyKVFinder/

Documentation for the Python API of the Schrodinger Maestro protein-editing tool suite can be found here:
https://learn.schrodinger.com/public/python_api/2025-1/overview.html
In addition, the two very basic example script files provided in this directory should help, as I have verified that they work and perform the correct function:
/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/example_script/mutate_m102g.py
- this script takes the sample structure, data/s0067_rec.pdb, and mutates one of its residues (Methionine 102) to a Glycine.
/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/example_script/change_m102_rotamer.py
- this script takes the sample structure, data/s0067_rec.pdb, and switches the rotamer state of three of its residues, Methionine 102, Leucine 118, and Valine 87.

# Relevant files and information

The experimental reference structure is:
/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/data/4W59_n-hexylbenzene.pdb

The correct ligand to select in the experimental reference structure has the name/three-letter code "3GZ". Please feel free to hardcode this selection in right now, but please provide the user the functionality to define which ligand three-letter code they'd like to select downstream.

The sample structure (single structure for now) we want to modify with this script is:
/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/data/s0067_rec.pdb

Please store modified structures in the /output/ directory, and ensure that the original sample structure name (e.g. "s0067_rec") is in the output filename.

Please store checkpoint files, structures or statistics, in the /checkpoints/ directory.

