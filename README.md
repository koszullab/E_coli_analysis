## Codes and functions for the multiscale analysis of the *Escherichia coli* 3C data ##

This page presents the different codes and functions developed to perform the analysis described in the article **Multiscale structuring of the *E. coli* chromosome by nucleoid-associated and condensin proteins** by Virginia S. Lioy, Axel Cournac, et al. The codes presented here should allow to reproduce the different graphs and figures from the main text and the supplementary data. 

### Table of contents

* [Dependencies](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building contacts map](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#building-contacts-map)
* [Correlation with other data](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#correlation-with-other-data)
* [Scalogram visualization tool](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#scalogram-vizulaisation-tool)
* [Directionality Index](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#directionality-index)
* [Correlation between transcription and contacts](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#correlation-between-transcription-and-contacts)
* [3D structure](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#3d-structure)
* [Ratio of contacts](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#ratio-of-contacts)


### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems, and necessitate:
#### *Python (>=2.7)*
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

FASTQ files of the reads were deposited in the NCBI database under the GEO accession number GSE107301. A SRA executable called fastq-dump from SRA can be used to extract and split both reads of pair-end sequences: 
```bash
fastq-dump library_identification --split-3 -O /path_to_a_directory
```


#### Alignment

We use the MG1655 reference genome (GenBank: U00096.2, total length 4639675), the program Bowtie2 and an iterative procedure (see for instance the one described in [DADE] (https://github.com/scovit/dade). We process the pairs of reads so only read with a mapping quality > 30 are retained. For instance:
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

# Filtering of pairs of reads that both have a Mapping Quality above 30
awk '{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}'  p1_p2_merged  > output_alignment_idpt.dat

# Removal of intermediar file
rm p1_p2_merged
```
At this stage, you have a file (output_alignment_idpt.dat.ind3) containing these information organized as such:
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

where chr1 corresponds to *Escherichi coli* genome. Each read is assigned to a restriction fragment as described in https://github.com/axelcournac/3C_tutorial (Cournac et al., MMiB, 2016).


## Building contacts map
To build the contact map and filter non-informative events, we use 3Cevents_MATRICE.py [`fragment_attribution.py`](python_codes/fragment_attribution.py):
```bash

python fragment_attribution.py /media/axel/RSG3/BACK_UP/axel/Bureau/python/fasta/ecoli/ HpaII output_alignment_idpt_BC76.dat

python library_events.py output_alignment_idpt_BC76.dat.indices BC76_MatP

python Matrice_Creator.py output_alignment_idpt_BC76.dat.indices.filtered 5000 BC76_MatP_37C
```
The first argument corresponds to the path of the file “output_alignment_idpt.dat.ind3”, which contains the alignment information. The second argument corresponds to the size of the bin (here: 5,000bp). The third argument is the name of the prefix for the file of the contacts maps. To plot and/or save the contact map, use [plot_mat_temp.py](python_codes/plot_mat_temp.py)
```bash
python plot_mat_temp.py mat_temp_WT_rep1_5000.txt WT_rep1_5000
```
![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/Ecolichromosomemap.png)
### Scalogram

The scalogram tool allows to vizualise the dispersion of the contacts signal along the spatial scales. The functions are implemented in the code [multi_scale_domainogram_FILES2_dom3_3plots.py](python_codes/multi_scale_domainogram_FILES2_dom3_3plots.py)

```bash
python multi_scale_domainogram_FILES2_dom3_3plots.py  mat_temp_WT_rep1_5000.txt WT_rep1_5000
```
![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/scalogram.png)

### Correlation with other data
To correlate 3C contacts and recombination previously generated in Valens et al., EMBO 2004, we used [recombination_3C_Correlation.py](python_codes/recombination_3C_Correlation.py)

![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/recombinationvs3C.png)

We use a similar code to compare the geometrical distances measured from microscopy and the ones extracted from 3D structure built with 3C contacts data:  [distance_structure.py](python_codes/distance_structure.py)

To correlate MSD (mean square displacements) from time lapse microscopy experiments and cumulative contacts signal, the code [correlation_MSD_3C.py](python_codes/correlation_MSD_3C.py) was used. 

### Directionality Index 

We computed at two different scales: at 400 kb scale (macrodomains) and 100 kb scale (CIDs). The computation was done on the correlation matrices (to attenuate fluctuations in the signals) as previously described (Dixon et al., 2012).

### Correlation between transcription and contacts

To put in evidence the correlation between transcription and 3C contacts at short range, we use the code  [correlation_transcription_3C.py](python_codes/correlation_transcription_3C.py)

We use a Gaussian function to smooth both signals. The following graphs are generated with this package:

![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/correlation_transcription_rna_olivier_3C_MQ30_3.png)

To plot the frequency of contact as a function of genomic distance for different groups of bins classified according to their transcription level, we first compute the distribution of log2 number of reads from the transcriptome for every bins. We divided the bins into 3 groups: poorly expressed bins with transcription level < 7, moderately expressed with transcription level > 7 and < 10 and highly expressed with transcription level > 10.
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
![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/histo_transcription_olivier.jpeg)

We then computed the genomic distance on the normalized contact maps (poor interacting bins were removed during this process) with:
```python 
import distance_law3
import scn

MAT_INT = loadtxt('mat_temp_WT_fused_4banks_5000.txt')
matscn1= scn.scn_func(MAT_INT,100);

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
This gives the following graph.

![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/dist_laws_transrciption_groups_rnaseq_olivier.jpeg)



### Ratio of contacts
We computed the ratio of contacts between mutants and the corresponding WT maps by averaging the contacts made at a certain distance for mutant and WT normalized contact maps and computing their log2 ratio. The computation is implemented in the code [RATIO_CONTACTS_2dtype.py](python_codes/RATIO_CONTACTS_2dtype.py)

```bash
python RATIO_CONTACTS_2dtype.py mat_temp_5000_WT.dat mat_temp_5000_MatP.dat MatP
```

![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/ratio_of_contacts.png)


To plot on short scales:

```python
s=s[range(s.shape[0]-1,-1,-1),:];
s2=s[range(0,8),:];
extent=(0, 928, 0, 80)
imshow( s2[range(s2.shape[0]-1,-1,-1),:],interpolation="none", extent=extent, vmin=-1.0, vmax=1.5,cmap="seismic")
colorbar();
```

Here an example for the short scales zoom of the ratio of contacts for HNS/WT:

![alt tag](https://github.com/koszullab/E.coli.analysis/blob/master/pictures/h-ns_vs_wt.png)


### Miscellaneous

To analyze plasmids contacts signal along the main Escherichia coli genome and plot the signal centered on matS sites, we used the R code [plasmid_genome_norm4.R](R_codes/plasmid_genome_norm4.R).

The enrichment analysis of the CIDs borders was carried out thanks to the R script [borders_statistics5.R](R_codes/borders_statistics5.R).


