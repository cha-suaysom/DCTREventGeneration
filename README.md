# DCTREventGeneration

1. Install Pythia8 http://home.thep.lu.se/~torbjorn/Pythia.html and Delphes (https://github.com/delphes/delphes)

2. Follow this instruction to get Delph+Pythia8 to work https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/Pythia8

./DelphesPythia8 examples/DCTREventGeneration/delphes_card_CMS.tcl examples/Pythia8/configNoLHE.cmnd output.root

root -b -x -q 'examples/DCTREventGeneration/mymass.C("output.root","output_det.txt","output_part.txt")'