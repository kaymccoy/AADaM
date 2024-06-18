### CITATION ###

To cite AADaM, please cite [link to be added].



### LICENSE ###

MIT License

Copyright (c) 2024 Katherine McCoy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



### INTRODUCTION ###

The Antibody-Antigen Dataset Maker (AADaM) is python script that takes a larger dataset downloaded from SAbDab, and uses it to create benchmark / testing datasets intended for ML-methods. It may also create complementary training datasets for ML methods that use antibody-antigen structures for training!

Structures accepted into the dataset will be required to have below a certain % CDR sequence identity with protein-binding antibody structures released past the training cutoff date (-d flag), to ensure each benchmark is fair with respect to testing ML methods, in that it will contain no close similarities to potential training data. Then within the benchmark, targets are disallowed from sharing a certain % or more sequence identity between heavy chain CDRs, light chain CDRs, and/or any antigen chains to remove redundancy from the benchmark (-cd flag), for the sake of removing representation bias. The and/or is determined by how strict you want to be (-cs flag).

The dataset produced may also be filtered my method (-m flag), resolution (-r flag), and non-natural residues (-nx flag). A whitelist of structures you want to have preference to be included may also be rpovided (-w flag).

The benchmark was filtered for various indications of quality: resolution, method used to determine the structure, gaps in the CDRs or antigen, length of the antigen, etc. (Fig S1, also see Materials & Methods).

Finally, there are two different ways to determine sequence ID %; either by using local alignment and requiring it to be 30 residues long or greater, or by using global alignment. 

When one structure "knocks out" another from consideration, the one with fewest missing residues within the CDRs and antigen chain(s) is given preference to remain within the dataset, and if both structures share the same number of missing residues (e.g. both have no missing residues), the structure with the shorter antigen sequence is selected. 



### SET UP ###

To set up AADaM, please first set up Mosaist, located on GitHub at Grigoryanlab/Mosaist.

Then provide your path to Mosaist's lib directory on the 7th line of AntibodyAntigenDatasetMaker.py, replacing "/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/lib"



### HELPER SCRIPTS ###

Also included in the repo are helper scripts to IMGT number and otherwise clean up antibody-antigen structures, to search antibody-antigen interfaces for structural motifs, to check the Neff of both paired and single chain MSAs, and to calculate interfacial pLDDT. These scripts may be useful for analyzing antibody-antigen models in your future projects! Those relying on Mosaist similarly require replacing the string for the Mosaist lib to work.
