# Amino_acid_composition_and_hydrophobicity_analysis

Project description:
This package is a simple set of bioinformatics to analyze the amino acid content and hydrophobicity 
profile of a given set of protein sequences in a fasta file. 

How to install the project:
This code only uses stock Matlab functions. This code should run on any operating system capable of 
running Matlab, but has been formally tested on Windows 10. No special hardware is required. If Matlab 
is already installed, installing the package should take less than 2 min. The code has been tested on 
Matlab R2020b.

Running the code:
The composition_hydrophobicity.m code takes a reference proteome fasta file, as well as a fasta file 
containing your test sequences as inputs. The code then calls on count_aa_in_string.m and 
count_aa_in_fasta.m functions to look at the amino acid composition in your test sequences as well as 
the reference proteome. Then, the code calculates the hydrophobicity profile of the test sequences. 
The results are displayed in plots and also stored in a structure called 
“composition_hydrophobicity_struct.”
An example of this type of analysis can be found in test_code.m. The code took about 18 seconds to run 
on a 13.4 KB genome on my machine.

Acknowledgments:
This work was supported by the joint efforts of The Michael J. Fox Foundation for Parkinson’s Research 
(MJFF) and the Aligning Science Across Parkinson’s (ASAP) initiative. MJFF administers the grants ASAP-
000282 and ASAP- 024268 on behalf of ASAP and itself.
