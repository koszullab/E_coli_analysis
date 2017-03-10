# Codes and functions for the multiscale analysis of the *Escherichia coli* 3C data 
This page presents the different codes and functions developped for the multiscale analysis of 3C data of the organism *Escherichia coli* presented in the article **Cooperation of a condensin complex and nucleoid-associated proteins for the multiscale conformation of a bacterial chromosome** by Virginia S. Lioy, Axel Cournac, Martial Marbouty, Stéphane Duigou, Julien Mozziconacci, Olivier Espéli, Frédéric Boccard, Romain Koszul.
The codes presented here should allow to exactly reproduce the various plots present in the main manuscript and supplementary data. 


For queries or help getting these running, you can send email or open an issue at the github repository.

### Table of contents

* [Dependencies](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#building-of-the-contacts-map)
* [Scalogram visualization tool](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#Scalogram vizulaisation tool)
* [Directionality Index tool tool](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#Directionality-Index-at-400-kb-scale-(macrodomains)-and-100-kb-scale-(CIDs))
* [Correlation between transcription and 3C](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#Correlation-between-transcription-and-3C-contacts)
* [3D structure](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#3D-structure)
* [Ratio of contacts](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#Ratio-of-contacts-between-mutant-and-corresponding-WT-along-the-spatial-scales)


### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `Pymol` / [Pymol](https://www.pymol.org/)
* `R` / [R](https://cran.r-project.org/)
* `Shrec3D` / [Shrec3D](https://sites.google.com/site/julienmozziconacci/)



## Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA can be used to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory

```bash
./fastq-dump SRR639031 --split-3 -O /run/media/axel/EColi_data/
```

#### Alignment

For the alignement step, we will use the sofware Bowtie2 and an iterative procedure like the one of [hiclib] (https://bitbucket.org/mirnylab/hiclib). 

We process the pairs of reads so that every read has a mapping quality superior to 30. 
Here, some lines that can be used to do this task:

```bash
#  Keeping only the columns of the sam file that contain necessary information:
awk '{print $1,$3,$4,$2,$5;}' p1.sam > p1.sam.0
awk '{print $1,$3,$4,$2,$5;}' p2.sam > p2.sam.0

# Sort according to the read identification to have both mates in the same order
# if sort does not have -V option try -d
sort -V -k1 p1.sam.0 > p1.sam.0.sorted
sort -V -k1 p2.sam.0 > p2.sam.0.sorted

# Pairing of both mates in a single file
paste p1.sam.0.sorted p2.sam.0.sorted > p1_p2_merged

# Removal of intermediar files
rm p1.sam.0.sorted
rm p2.sam.0.sorted

# Filtering of paires of reads that both have a Mapping Quality above 30
awk '{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}'  p1_p2_merged  > output_alignment_idpt.dat

# Removal of intermediar file
rm p1_p2_merged
```
At this stage, you should have a file containing these information:
```
chr1 104180 16 chr1 104057 0
chr1 3570510 16 chr1 3570450 0
chr1 4255981 0 chr1 4256104 16
chr1 159457 16 chr1 159370 0
chr1 4113710 16 chr1 4113584 0
chr1 4259849 16 chr1 4259818 0
chr1 3753874 0 chr1 3754001 16
chr1 2856270 16 chr1 2856124 0
chr1 4134782 16 chr1 4134678 0
```

chr1 corresponds here to the chromosome of *Escherichi coli* genome. We used  MG1655 reference genome (GenBank: U00096.2, total length 4639675). We name this file output_alignment_idpt.dat.ind3.

We then assigned eahc read to its corresponding restriction fragment as described previously in [https://github.com/axelcournac/3C_tutorial](https://github.com/axelcournac/3C_tutorial). 



## Building of the contacts map
To build the contact map and and filtered the non informative events, we use the python code 3Cevents_MATRICE.py [`3Cevents_MATRICE.py`](python_codes/3Cevents_MATRICE.py):
```bash
python 3Cevents_MATRICE.py output_alignment_idpt.dat.ind3 5000 WT_rep1_5kb
```
The first argument corresponds to the path of the file named output_alignment_idpt.dat.ind3 containing the informations of the mapped pairs of reads. 
The second argument is the size of the bin, here 5000 bp. We generally use this resolution for the whole study which is a good compromise between resolution and signa robustness.
The third argument is the name of the prefixe for the file of the contacts maps. 

The generated picture is and corresponds to the first figure of the manuscript. 
To plot the contact map and save in various formats, we use the python code [plot_mat_temp.py](python_codes/plot_mat_temp.py)
```bash
python plot_mat_temp.py mat_temp_WT_rep1_5000.txt WT_rep1_5000
```

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/WT_rep1_5000_Natural.png)


## Directionality Index at 400 kb scale (macrodomains) and 100 kb scale (CIDs)

We computed and plotted the Directional Index tool as previously described. 

## Scalogram vizulaisation tool 

The scalogram tool allows to vizualise the dispersion of the contacts signal along the spatial scales. The functions are implemented in the code [multi_scale_domainogram_FILES2_dom3_3plots.py](python_codes/multi_scale_domainogram_FILES2_dom3_3plots.py)

```bash
python multi_scale_domainogram_FILES2_dom3_3plots.py  mat_temp_WT_rep1_5000.txt WT_rep1_5000
```
![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/WT_rep1_5000_DOM.jpeg)

The second argument is the prefixe for names of output files.

## Correlation between transcription and 3C contacts
To put in evidence the correlation bewteen transcription and 3C contacts at short range, we use the python code [correlation_transcription_3C.py](python_codes/correlation_transcription_3C.py)

We use gaussian function to filter both signals and attenuate the fluctuations. 
This code generates the following figure:

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/correlation_transcription_rna_olivier_3C.png)

## 3D structure
For the construction of 3D structure, we processed the matrice by removing the outliers elements. We computed the genomic distance law and removed points outside the mean + 2 std using the function 'filter_dist'.
We use the algorithm Shrec3D with the modification that the lwa to convert contact frequencies into geometricla distance is d=(1/f^0.5).
We use home made pymol script to generate the 3D picture including the strong contact of the 3C library. 

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/EColi_genome_local_constrain252.png)


## Ratio of contacts between mutant and corresponding WT along the spatial scales 
We computed the ratio of contacts between mutant and correponding WT by averaging the concacts made at a certain distance for mutant and WT normalised contact maps and then we took the log2 ratio. 
The computation is implemented in the code [RATIO_CONTACTS_2dtype.py](python_codes/RATIO_CONTACTS_2dtype.py)

```bash
python RATIO_CONTACTS_2dtype.py mat_temp_5000_WT.dat mat_temp_5000_MatP.dat MatP
```

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/MatP_scalogram_RATIO_2D_GAUSSIAN_2_LOG_SEISMIC.png)




