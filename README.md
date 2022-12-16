# ROS-ALE
Generates all figures for the manuscript "Lab evolution, transcriptomics, and modeling reveal mechanisms of paraquat tolerance" by Kevin Rychel, et al.

## Additional Repositories
1. [Modulome-Workflow](https://github.com/avsastry/modulome-workflow): Aligns RNA sequencing reads to the genome and computes iModulons. input/precise1k_kr.json was generated using this repository.
2. [PyModulon](https://github.com/SBRG/pymodulon): Class and functions for analyzing iModulons. Functions for reading input/precise1k_kr.json, making iModulon-related plots, and generating iModulonDB.org pages are in this repository.
3. [OxidizeME](https://github.com/laurenceyang33/oxidizeme): This is the Metabolic and Expression (ME) model that incorporates ROS stress. It is used to generate files in input/ME/, which are used to generate parts of Figure 5.

## Large Input File Links
Some files are too large to share on GitHub. If you would like to run this code, please clone this repository and add the following files from Dropbox in their respective locations:
1. [input/precise1k_kr.json](https://www.dropbox.com/s/6yxzdz0odjnxrtj/precise1k_kr.json?dl=0): Main iModulon object, containing transcriptomic data.
2. [input/pitA/pitA_data.json](https://www.dropbox.com/s/upw3vcf4rpqq29v/pitA_data.json?dl=0): PitA mutant strains iModulon object, which is used in Figure 4
3. [input/genome_coverage/1_0_reference.gff](https://www.dropbox.com/s/amz34udjp44kthv/1_0_reference.gff?dl=0): GFF for genome coverage in the 1_0 strain, used in Figure 3.
4. [input/genome_coverage/3_0_reference.gff](https://www.dropbox.com/s/557myw8vozjdu2j/3_0_reference.gff?dl=0): GFF for genome coverage in the 3_0 strain, used in Figure 3.
