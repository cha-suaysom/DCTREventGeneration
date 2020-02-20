# DCTREventGeneration

1. Install Pythia8 http://home.thep.lu.se/~torbjorn/Pythia.html and Delphes (https://github.com/delphes/delphes)

2. Follow this instruction to get Delph+Pythia8 to work https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/Pythia8

3. Put this folder in the same directory that Delphes is installed

4. Run the following commands

	1. `./DelphesPythia8 DCTREventGeneration/delphes_card_CMS.tcl DCTREventGeneration/configNoLHE.cmnd output.root`

	2. `root -b -x -q 'DCTREventGeneration/mymass.C("output.root","output_det.csv","output_part.csv","output_hadron.csv")'`	