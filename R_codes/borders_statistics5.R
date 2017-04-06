# script to assign several feature to a set of borders based on DI analysis
# 5: october 2016 to integrate other GO groups 

`%ni%` = Negate(`%in%`)

# input file for the borders: 
borders=read.table("/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/enrichment_borders/bins_borders_files/borders_detected_WT_BC70_5kb.txt")
#borders=read.table("/home/axel/Bureau/enrichment_borders/bins_borders_files/borders_detected_WT_MM30C_5kb.txt");
#borders=read.table("/home/axel/Bureau/enrichment_borders/borders_files/borders_detected_WT_LB30C_5kb.txt");
#borders=read.table("/home/axel/Bureau/enrichment_borders/borders_files/borders_detected_WT_LB37C_5kb.txt");

#borders=read.table("/home/axel/Bureau/DI_all_banks_100kb/replicats_MM_30C/borders_detected_MM_30C_NextSeq_dec2014.txt");
#borders=read.table("/home/axel/Bureau/DI_phase_stat/borders_detected_GTCT_Coli_Phae_stat_rep1.txt");
#borders=read.table("/home/axel/Bureau/enrichment_borders/bins_borders_files/borders_bins5kb_DI100kb_Stationary.txt")

#FILE="test";
#FILE="MM_30C_NextSeq_dec2014";
#FILE="DI100kb_Stationary"
FILE="DI100kb_BC70_5kb"

nb_borders = dim(borders)[1]
borders=borders$V1

# Biological features:
gene=read.table("/home/axel/Bureau/bacteries_project/coli/genes_K12_UCSC.dat2")

gene_expr = read.table("/home/axel/Bureau/bacteries_project/coli/genes_expression_MM.txt2")
#gene_expr = read.table("/home/axel/Bureau/bacteries_project/coli/genes_expression_LB.txt2")

on_genes=subset(gene_expr,gene_expr$V2 >250);
gene_sorted = gene_expr[with(gene_expr, order(-gene_expr$V2) ), ]
pc=10./100;
higly_genes = gene_sorted[ 1:floor(dim(gene_sorted)[1] * pc),];  # we take the 10% most expressed genes from Oliver data
rna_genes= read.table("/home/axel/Bureau/bacteries_project/coli/rrna.txt2");
higly_genes=rbind(higly_genes,rna_genes);

gene_sorted = gene_expr[with(gene_expr, order(gene_expr$V2) ), ];
pc=10./100;
poorly_genes = gene_sorted[ 1:floor(dim(gene_sorted)[1] * pc),];
# to remove rrna genes that can be in poorly expressed:
poorly_genes = subset(poorly_genes,  poorly_genes$V1 %ni% higly_genes$V1)

is = read.table("/home/axel/Bureau/bacteries_project/coli/me/IS_ALL.txt")

heEPOD=read.table("/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/Protein-occupancy/EPOD_enrichment/EPOD-technique/heEPOD.txt")
tsEPOD=read.table("/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/Protein-occupancy/EPOD_enrichment/EPOD-technique/tsEPOD.txt")

is = tsEPOD

SRP = read.table("/home/axel/Bureau/bacteries_project/coli/genes_SRP.txt");
SRP_on = subset(SRP, SRP$V1 %in% on_genes$V1);
#SRP_on = subset(SRP, SRP$V1 %in% on_genes$V1 & SRP$V1 %ni% higly_genes$V1);


SECB=read.table("/home/axel/Bureau/bacteries_project/coli/genes_secB.txt");
SECB_on = subset(SECB, SECB$V1 %in% on_genes$V1);

HTG=read.table("/home/axel/Bureau/bacteries_project/coli/me/ecoli.HGT.txt3");
HTG_on = subset(HTG, HTG$V1 %in% on_genes$V1);

#SECB_on = HTG


#  Computation:
computation(borders)

# # random realisations: 
# max=0;
# vect1=vector();vect2=vector();vect3=vector();
# for(r in 0:100)
# {
# object =computation(floor(runif(nb_borders, 0, 928) ) )
# vect1 = c(vect1,object[2]); vect2 = c(vect2,object[3]); vect3 = c(vect3,object[4]); 
# }

#  function to compute and write into file:
computation <- function(pos)   # take a vector of positions to compute different number of highly, IS, SRP 
{
  nb_pos=length(pos)
  BIN = 5000;
  nb_high=0;
  nb_is=0;
  nb_srp=0;
  nb_poor=0;
  nb_secb=0;

  df = data.frame("Bin",";","Pos1",";","Pos2",";","genes",";","10% highly expr.",";","IS",";","SRP",";","10 poorly expr.",";","SecB");
  write.table(df, file = FILE, append = TRUE,quote = F,row.names = F,col.names = F);
  for(i in 1:nb_pos )
  {
    p1=(pos[i]-1)*BIN
    p2=(pos[i]+2)*BIN
    g = subset(gene$V1,gene$V2> p1 & gene$V3 < p2);
    g2="";
    for(j in g)  {g2=paste(g2,j,sep = " ") }
    
    hg = subset(g,g %in% higly_genes$V1);
    hg2=" ";
    for(j in hg)  {hg2=paste(hg2,j,sep = " ") }
    if( length(hg) > 0) {nb_high= nb_high+1};
    
    pg = subset(g,g %in% poorly_genes$V1);
    pg2=" ";
    for(j in pg)  {pg2=paste(pg2,j,sep = " ") }
    if( length(pg) > 0) {nb_poor= nb_poor+1};
    
    i1 = subset(is$V1, is$V1 > p1 & is$V1 < p2);
    if( length(i1) > 0) {i2="yes";nb_is=nb_is+1;} else {i2="no";} 
    
    sg = subset(g,g %in% SRP_on$V1);
    sg2=" ";
    for(j in sg)  {sg2=paste(sg2,j,sep = " ") }
    if( length(sg) > 0) {nb_srp=nb_srp+1};
    
    sec_g = subset(g,g %in% SECB_on$V1);
    sec_g2=" ";
    for(j in sec_g)  {sec_g2=paste(sec_g2,j,sep = " ") }
    if( length(sec_g) > 0) {nb_secb=nb_secb+1};
    
    df <- data.frame( borders[i],";", p1,";", p2,";", g2,";",hg2,";",i2,";",sg2,";",pg2 ,";",sec_g2);
    write.table(df, file = FILE, append = TRUE,quote = F,row.names = F,col.names = F);
  }
  object= c(nb_pos, nb_high, nb_is, nb_srp, nb_poor,nb_secb) ;
  return(object)
}

#  Counting the number of occurences for biological features:
ob=computation(borders)   #  to have the number of occurences for the determined borders
oe=computation(1:928)     #  to have expectations. 

#  Fisher tests:
mat <- matrix( c(ob[2], nb_borders-ob[2], oe[2] ,928-oe[2]), nrow = 2);  tes <- fisher.test(mat);tes;p1=tes$p.value;
mat <- matrix( c(ob[3], nb_borders-ob[3], oe[3] ,928-oe[3]), nrow = 2);  tes <- fisher.test(mat);tes;p2=tes$p.value;
mat <- matrix( c(ob[4], nb_borders-ob[4], oe[4] ,928-oe[4]), nrow = 2);  tes <- fisher.test(mat);tes;p3=tes$p.value;
mat <- matrix( c(ob[5], nb_borders-ob[5], oe[5] ,928-oe[5]), nrow = 2);  tes <- fisher.test(mat);tes;p4=tes$p.value;
mat <- matrix( c(ob[6], nb_borders-ob[6], oe[6] ,928-oe[6]), nrow = 2);  tes <- fisher.test(mat);tes;p5=tes$p.value;
#


# Computation around the borders:
v1=vector();v2=vector();v3=vector();v4=vector();v5=vector();
for(b in -20:20)
{
  print(b); 
  object = computation(borders+b);
  v1 = c(v1,object[2]); v2 = c(v2,object[3]); v3 = c(v3,object[4]);v4 = c(v4,object[5]); v5 = c(v5,object[6]);  
}  




# Computation around the borders:
v1=vector()
v2=vector()
v3=vector()
v4=vector()
v5=vector();
for(b in -20:20)
{
print(b); 
object = computation_GO(borders+b);
v1 = c(v1,object[2]); v2 = c(v2,object[3]); v3 = c(v3,object[4]);v4 = c(v4,object[5]); v5 = c(v5,object[6])
}  


