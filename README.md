# Codes and functions for the multiscale analysis of the *Escherichia coli* 3C data 
This page presents the different codes and functions developped for the multiscale analysis of 3C data of the organism *Escherichia coli* presented in the article **Cooperation of a condensin complex and nucleoid-associated proteins for the multiscale conformation of a bacterial chromosome** by Virginia S. Lioy, Axel Cournac, Martial Marbouty, Stéphane Duigou, Julien Mozziconacci, Olivier Espéli, Frédéric Boccard, Romain Koszul.
The codes presented here should allow to exactly reproduce the various plots present in the main manuscript and supplementary data. 


For queries or help getting these running, you can send email or open an issue at the github repository.

### Table of contents

* [Dependencies](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#building-of-the-contacts-map)
* [Correlation with other data](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#correlation-with-other-data)
* [Scalogram visualization tool](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#scalogram-vizulaisation-tool)
* [Directionality Index](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#directionality-index)
* [Correlation between transcription and contacts](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#correlation-between-transcription-and-contacts)
* [3D structure](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#3d-structure)
* [Ratio of contacts](https://github.com/axelcournac/EColi_analysis/blob/master/README.md#ratio-of-contacts)


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
### Scalogram vizulaisation tool

The scalogram tool allows to vizualise the dispersion of the contacts signal along the spatial scales. The functions are implemented in the code [multi_scale_domainogram_FILES2_dom3_3plots.py](python_codes/multi_scale_domainogram_FILES2_dom3_3plots.py)

```bash
python multi_scale_domainogram_FILES2_dom3_3plots.py  mat_temp_WT_rep1_5000.txt WT_rep1_5000
```
![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/WT_rep1_5000_DOM.jpeg)

The second argument is the prefixe for names of output files.

### Correlation with other data
To correlate 3C contacts and recombination prevuously generated in Valens et al., EMBO 2004, we used the python code [recombination_3C_Correlation.py](python_codes/recombination_3C_Correlation.py)

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/recombination_vs_3C.png)

We use a similar code to confront geometrical distances mesaured from microscopy and the ones extracted from 3D structure built with 3C contacts data: [distance_structure.py](python_codes/distance_structure.py)

To correlate MSD (Mean square Displacements) from time lapse microscopy experiments and cumulative contacts signal, the code [correlation_MSD_compaction_Marco.py](python_codes/correlation_MSD_compaction_Marco.py) was used. 

### Directionality Index 

We computed at two different scales: at 400 kb scale (macrodomains) and 100 kb scale (CIDs). The computation was done on the correlation matrices (to attenutate fluctuations in the signals) as previously described. 

### Correlation between transcription and contacts

To put in evidence the correlation bewteen transcription and 3C contacts at short range, we use the python code [correlation_transcription_3C.py](python_codes/correlation_transcription_3C.py)

We use gaussian function to smooth both signals and attenuate the fluctuations. 
This code generates the following figure:

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/correlation_transcription_rna_olivier_3C_MQ30_3.png)

To plot the genomic distance law for different groups of bins classified according to their transcription level, we first compute the distribution of log2 number of reads from the transcriptome for every bins. We divided the bins into 3 groups: poorly expressed bins with transcription level < 7, moderately expressed with transcription level > 7 and < 10 and highly expressed with transcription level > 10. 

```python 
chip=loadtxt("hist_EV-4_TGACCA_L002_R1_001.MQ30.hist5000")
#plot(log(chip[:,1]) );

c,h,w=my_histo.my_histo( log(chip[:,1]), 100);
plt.bar(c, h, align='center', width=w)
plt.bar(c[c<7], h[c<7], align='center', width=w,color="green");
plt.bar(c[c>10], h[c>10], align='center', width=w,color="red");

plt.xlabel("Transcription level (log2 of number of reads)");
plt.ylabel("Number of occurences (bins of 5kb))");
plt.show();
```
![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/histo_transcription_olivier.jpeg)

We then computed the genomic distance on the normalised matrice (poor interacting bins were removed during this process) with the following lines of code:
```python 
import distance_law3
import scn

MAT_INT = loadtxt('mat_temp_WT_fused_4banks_5000.txt')
matscn1= scn.scn_func(MAT_INT,15000);

# 2)  with rna-seq from Olivier Espeli 
# very expressed:
indices_highly=np.where( log(chip[:,1]) > 10);
#indices_highly=indices_highly[0];
mat_highly= MAT_INT[indices_highly[0],:];
mat_highly.shape;

dih=distance_law3.dist_law(mat_highly,indices_highly[0]);

# poorly expressed:
indices_poor=np.where( log(chip[:,1]) < 7);
#indices_poor=indices_poor[0];
mat_poor= MAT_INT[indices_poor[0],:];
mat_poor.shape
dip=distance_law3.dist_law(mat_poor,indices_poor[0]);

# normaly expressed:
#indices_norm=np.where( (log(chip) > 6) and (log(chip) < 9) );
indices_norm = np.where( logical_and( (log(chip[:,1]) > 7) , (log(chip[:,1]) < 10) ) );

#indices_poor=indices_poor[0];
mat_norm= MAT_INT[indices_norm[0],:];
mat_norm.shape
din=distance_law3.dist_law(mat_norm,indices_norm[0]);

plot(dih,color="red",label="very transcribed",linewidth=4.0);
plot(din,color="blue",label="moderately transcribed",linewidth=4.0);
plot(dip,color="green",label="poorly transcribed",linewidth=4.0);
xlabel("Genomic distance (in bins of 5kb)");
ylabel("Number of reads per possible distance")
legend();
grid();
```
This gives the following graph showing that transcription process increases contacts at short scales.

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/dist_laws_transrciption_groups_rnaseq_olivier.jpeg)



### 3D structure
For the construction of 3D structure, we processed the matrice by removing the outliers elements. We computed the genomic distance law and removed points outside the mean + 2 std using the function 'filter_dist'.
We use the algorithm Shrec3D with the modification that the lwa to convert contact frequencies into geometricla distance is d=(1/f^0.5).
We use home made pymol script to generate the 3D picture including the strong contacts by adding links between monomers. 
```python
bg_color white

#cmd.load("/home/axel/Bureau/script_pymol/3Dcoor_ecoli_5kb_3_2.xyz","mov")
#cmd.load("/home/axel/Bureau/files_for_shrec3D/3Dcoor_ecoli_3.xyz","mov")

cmd.load("/home/axel/Bureau/files_for_shrec3D/3Dcoor_ecoli_WT_Hiseq_WT_MM_30C_BC164_2.xyz","mov")
cmd.load("/home/axel/Bureau/files_for_shrec3D/3Dcoor_ecoli_matp_2.xyz","mov")


#n = 864
n = 908

set sphere_scale, 2, (all)
show sphere
for idx in range(0,n):cmd.bond('index %s' % str(idx), 'index %s' % str(idx+1))
# for idx in range(0,n):cmd.unbond('index %s' % str(idx), 'index %s' % str(idx+1))

#inFile = open("/home/axel/Bureau/script_pymol/cumulative3C_200kb.txt", 'r')
inFile = open("/home/axel/Bureau/files_for_shrec3D/cumulative_signal_200kb_matrice_WT1_for_shreck.txt")

inFile = open("/home/axel/Bureau/files_for_shrec3D/cumulative_signal_200kb_WT_MM_30C_BC164_SCN_despeckled.txt")
inFile = open("/home/axel/Bureau/files_for_shrec3D/cumulative_signal_200kb_matP.txt01")

# create the global, stored array
stored.newB = []
# read the new B factors from file
for line in inFile.readlines(): stored.newB.append( float(line) )
# close the input file
inFile.close()
# clear out the old B Factors
alter mov, b=0.0
# update the B Factors with new properties
alter mov and n. A, b=stored.newB.pop(0)
# color the protein based on the new B Factors of the alpha carbons
cmd.spectrum("b", "blue_white_red", "mov and n. A")

set sphere_scale, 3, (all)
show sphere
set_bond stick_radius, 1.0, all
```

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/EColi_genome_local_constrain252.png)


### Ratio of contacts
We computed the ratio of contacts between mutant and correponding WT by averaging the concacts made at a certain distance for mutant and WT normalised contact maps and then we took the log2 ratio. 
The computation is implemented in the code [RATIO_CONTACTS_2dtype.py](python_codes/RATIO_CONTACTS_2dtype.py)

```bash
python RATIO_CONTACTS_2dtype.py mat_temp_5000_WT.dat mat_temp_5000_MatP.dat MatP
```

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/MatP_scalogram_RATIO_2D_GAUSSIAN_2_LOG_SEISMIC.png)


To plot on short scales:

```python
s=s[range(s.shape[0]-1,-1,-1),:];
s2=s[range(0,8),:];
extent=(0, 928, 0, 80)
imshow( s2[range(s2.shape[0]-1,-1,-1),:],interpolation="none", extent=extent, vmin=-1.0, vmax=1.5,cmap="seismic")
colorbar();
```

Here an example for the short scales zoom of the ratio of contacts for HNS/WT:

![alt tag](https://github.com/axelcournac/EColi_analysis/blob/master/pictures/HNS_zoom.png)

