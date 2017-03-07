# Codes and functions for the multiscale analysis of the Escherichia coli 3C data 
This page presents and explains the different codes and functions developped for the analysis of 3C data of the organism *Escherichia coli* presented in the article **Cooperation of a condensin complex and nucleoid-associated proteins for the multiscale conformation of a bacterial chromosome** by Virginia S. Lioy, Axel Cournac, Martial Marbouty, Stéphane Duigou, Julien Mozziconacci, Olivier Espéli, Frédéric Boccard, Romain Koszul.
The codes presented here allow to exactly reproduce the various plots present in the main manuscript and supplementary data. 
For queries or help getting these running, you can send email or open an issue at the github repository.

### Table of contents

* [Dependencies](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#building-of-the-contacts-map)
* [Normalization of the data](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#normalization-of-the-data)
* [Computation of genomic distance law](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-genomic-distance-law)
* [Computation of correlation matrices](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-correlation-matrices)
* [Directional Index tool to detect TADs](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#directional-index-tool-to-detect-tads)
* [Decomposition into eigen vectors](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#decomposition-into-eigen-vectors)
* [Use of sparse formalism](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#use-of-sparse-formalism)
* [Miscellaneous](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#miscellaneous)

### Dependencies

Scripts will run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `SRA tool` / [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)
* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
