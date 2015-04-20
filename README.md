# Assessing the aligner's behavior 

This repository contains a number of Python scripts to assess the behavior of a sequence aligner and is made public for our participation in the ICGC-TCGA DREAM Genomic Mutation Calling Challenge (syn312572). 

***

## Dependencies

The project depends for Python on the following packages: 

* _PySAM_ (see https://code.google.com/p/pysam/) for working with BAM/SAM files and
* _PyVCF_ (see https://github.com/jamescasbon/PyVCF) for working with VCF files. 

PySAM requires the installation of 

* _SAMtools_ (see http://samtools.sourceforge.net)

In addition, SimSeq should be installed: 

* _SimSeq_ (see https://github.com/jstjohn/SimSeq)

*** 

## Directory structure 

The repository consists of the following directories: 

* `bin/` - contains the Python scripts. Each script is equipped with a usage description. 

* `python/` - contains Python code that is re-used throughout some of the scripts in the `bin`-folder.

* `test/` - contains a small randomly generated test data set. 

***

## Contact

Louis Dijkstra

__E-mail__: louisdijkstra (at) gmail.com