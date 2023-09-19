[![DOI](https://zenodo.org/badge/293524912.svg)](https://zenodo.org/badge/latestdoi/293524912)
[![DOI:10.7717/peerj.16028](http://img.shields.io/badge/DOI-10.7717/peerj.16028-02a7fc.svg)](https://doi.org/10.7717/peerj.16028)

# **PHFinder: Assisted detection of point heteroplasmy in Sanger sequencing chromatograms.**

## Overview

**PHFinder** is an open-source Python tool designed to automatically assist on the detection of point heteroplasmy in Sanger sequencing chromatograms (AB1 files). This tool was designed with the aim of reducing the amount of files to be manually checked, especially on large datasets and set a less subjective definition of what double peaks can be considered an heteroplasmy.

A more detail description of PHFinder and its functioning can be found in the published paper:



## Installation

PHFinder consists of a main Python script and a Bash script. It requires to function the following dependencies:

+ [**Biopython**](https://biopython.org/)
+ [**Bowtie2**](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

After installation of these dependencies, to run PHFinder it's only necessary to copy the scripts in the desired working directory (and if necessary, make the files executable command: `chmod u+x PHFinder.py Bash_h.sh`).

It is also possible to run PHFinder in any directory without copying the scripts by downloading the folder PHFinder_Install and running inside the command `sudo python setup.py install`.  After this is done PHFinder can be used in any directory with the command `PHFinder.py`.

PHFinder scripts can be run on UNIX systems.


## Instructions for Analysis

PHFinder has been optimized to be run in a directory where the chromatograms are divided in several subfolders (see [**ExampleDataDirectory**](https://github.com/MSuarezMenendez/PHFinder/tree/master/ExampleDataDirectory), data also available [**here**](https://github.com/MSuarezMenendez/PHFinder/tree/master/ExampleDataDirectory)). But it will analyse any chromatogram present in the working directory (no matter the directory structure).

The arguments to use PHFinder in the terminal are the following:

+ **-o** Name of the output directory
+ **-f** Indicates that FASTQ files are to be extracted from the AB1 files
+ **-a** Fasta file to be used for the alignment
+ **-d** Enables the detection of heteroplasmy in the AB1 files
+ **-r** Ratio threshold for the detection of heteroplasmy (Default: 30)
+ **-q** Average quality threshold for the analysed region (Default: 40)
+ **-s** Analysis starting position (based on provided reference)
+ **-e** Analysis ending position (based on provided reference)
+ **-r** Secondary ratio thresholds for the detection of heteroplasmy (Default: 0.4) Can be deactivated by setting a high value (e.g., 100)
+ **-t** Test run **\***

Example run:

`PHFinder.py -o "Outputfolder" -f -a "Reference.fasta" -d -r 15 -q 40 -s 25 -e 350`

This command will extract fast files for every AB1 file present in the working directory, align them to the provided reference and look for double peaks (small peak at least 15% the size of the main one), from position 25 of the reference to position 350 if the average quality on this region is at least 40.

**\*** The test run option will automatically run 64 different combination of threshold indexes. To run this option a file named [**List_HP.csv**](https://github.com/MSuarezMenendez/PHFinder/tree/master/ExampleDataDirectory/List_HP.csv) with already known heteroplasmies in the dataset has to be present in the working directory. This allows to test what set of thresholds would detect most of the known heteroplasmies with least amount of ab1 files to go through.

## Citation

Suárez Menéndez M, Rivera-León VE, Robbins J, Berube M, Palsbøll PJ. 2023. PHFinder: assisted detection of point heteroplasmy in Sanger sequencing chromatograms. PeerJ 11:e16028 https://doi.org/10.7717/peerj.16028

## Version

The current release of **PHFinder** is [**v1.0**].


## License

GNU General Public License v3.0

Copyright (C) 2022 Marcos Suarez

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (version 3 of the License)

This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU General Public License for more details.

## Contact

PHFinder is written by Marcos Suarez (m.suarez.menendez@rug.nl)
