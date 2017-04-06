#  Script to analyse pattern of interactions between plasmid (+-MatS) and main genome
#  norm: to normalise with replication signal 

require("matrixStats")

# Miseq data:
#d=read.table("output_alignment_idpt.dat.ind3");
#d=read.table("/home/axel/Bureau/data/vick/exp4_dec2014_MiSeq/BC164_GTGT_MG1655_pGM2-5xmatS_30C_MM_REDONE/tmp/output_alignment_idpt.dat.ind3");

# Nextseq data:
#d11=read.table("/home/axel/Bureau/data/vick/exp7_juin2015_NextSeq/seqs/BC162_MG1655_pGM2_30C_MM/tmp/output_alignment_idpt.dat.ind3");
#d=read.table("/home/axel/Bureau/data/vick/exp7_juin2015_NextSeq/seqs/BC164_MG1655_pGM2-5xmatS_30C_MM/tmp/output_alignment_idpt.dat.ind3");

# Fusion of the two NextSeqs:
#d1=read.table("/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC162_GACT_MG1655_pGM2_30C_MM/tmp/output_alignment_idpt.dat.ind3");
#d2=read.table("/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC164_GTGT_MG1655_pGM2-5xmatS_30C_MM/tmp/output_alignment_idpt.dat.ind3");

# Miseq dec2014:
#d1=read.table("/home/axel/Bureau/data/vick/exp4_dec2014_MiSeq/BC162_GACT_MG1655_pGM2_30C_MM_REDONE/tmp/output_alignment_idpt.dat.ind3");
#d2=read.table("/home/axel/Bureau/data/vick/exp4_dec2014_MiSeq/BC164_GTGT_MG1655_pGM2-5xmatS_30C_MM_REDONE/tmp/output_alignment_idpt.dat.ind3");

# Miseq Feb 2016:
#d1=read.table("/home/axel/Bureau/data/vick/exp8_february2016_NextSeq/seqs/BC162_GACT_zapB_pGM2/tmp/output_alignment_idpt.dat.ind3");
#d2=read.table("/home/axel/Bureau/data/vick/exp8_february2016_NextSeq/seqs/BC110_CACT_zapB_pGM2-matS/tmp/output_alignment_idpt.dat.ind3");

#d1=read.table("/home/axel/Bureau/data/vick/exp8_february2016_NextSeq/seqs/BC164_GTGT_matPC20_pGM2/tmp/output_alignment_idpt.dat.ind3");
#d2=read.table("/home/axel/Bureau/data/vick/exp8_february2016_NextSeq/seqs/BC172_CGGT_matPC20_pGM2-matS/tmp/output_alignment_idpt.dat.ind3");

# NextSeq March 2016:
#d1=read.table("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/Ecoli_NextSeq_march_2016/BC164_GTGT_matPC20_pGM2/tmp/output_alignment_idpt.dat.ind3");
#d2=read.table("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/Ecoli_NextSeq_march_2016/BC172_CGGT_matPC20_pGM2-matS/tmp/output_alignment_idpt.dat.ind3");

#d1=read.table("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd1/Ecoli_NextSeq_march_2016/BC162_GACT_zapB_pGM2/tmp/output_alignment_idpt.dat.ind3");
#d2=read.table("/home/axel/Bureau/data/vick/exp10_April2016_MiSeq/116-145_S0_all_R1_2_001.fastq.pcrfree.out.MQ40.23");

# Fusion of the two NextSeqs (december 2016 filtered)
d1=read.table("/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC162_GACT_MG1655_pGM2_30C_MM/tmp/output_alignment_idpt.dat.ind3.fused2.dat3.indices.filtered2");
d2=read.table("/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC164_GTGT_MG1655_pGM2-5xmatS_30C_MM/tmp/output_alignment_idpt.dat.ind3.fused2.dat3.indices.filtered2");


#------------------------------------------------------------------
par(mfrow=c(1,1), cex = 0.5)
dev.off();
BIN = 1000
LEN_GEN = 4639675
#--------------------------------------------------------------

# for general coverage
dpg1 = subset(d1, d1$V1 == 1)
dpg2 = subset(d1, d1$V5 == 1)
dpg_fusioned = c(dpg1$V2,dpg2$V6);
h1 = hist(dpg_fusioned,breaks = seq(1,LEN_GEN+BIN,by = BIN), 
          xlab = "Position along the main genome (in bp)",
          ylab = "Number of contacts with plasmid pGM2",
          main =  paste("Contacts between plasmid and main genome with bin of",BIN,"bp") ,
          col = "deepskyblue4", 
          #      xlim =c(1e+06,2e+06) );
          cex = 0.4
)

# for general coverage
dpg1 = subset(d2, d2$V1 == 1 )
dpg2 = subset(d2, d2$V5 == 1 )
dpg_fusioned = c(dpg1$V2,dpg2$V6)

h2 = hist(dpg_fusioned,breaks = seq(1,LEN_GEN+BIN,by = BIN),
          xlab = "Position along the main genome (in bp)",
          ylab = "Number of contacts with plasmid pGM2",
          main =  paste("Contacts between plasmid and main genome with bin of",BIN,"bp") ,
          col = "deepskyblue4", 
          #      xlim =c(1e+06,2e+06)
)

cover1 = h1$counts
cover2 = h2$counts

#---------------------------------------------------------------------------------------------------
dpg = subset(d1, d1$V1 != d1$V5 );  # for the interactions with the plasmid
dpg1 = subset(dpg, dpg$V1 == 1);
dpg2 = subset(dpg, dpg$V5 == 1);
dpg_fusioned = c(dpg1$V2,dpg2$V6);
h1 = hist(dpg_fusioned, breaks = seq(1,LEN_GEN+BIN,by = BIN), 
          xlab = "Position along the main genome (in bp)",
          ylab = "Number of contacts with plasmid pGM2",
          main =  paste("Contacts between plasmid and main genome with bin of",BIN,"bp") ,
          col = "deepskyblue4", 
          #      xlim =c(1e+06,2e+06) );
          cex = 0.5
)

dpg = subset(d2, d2$V1 != d2$V5 );   # for the interactions with the plasmid
dpg1 = subset(dpg, dpg$V1 == 1 );
dpg2 = subset(dpg, dpg$V5 == 1 );
dpg_fusioned = c(dpg1$V2,dpg2$V6);

h2 = hist(dpg_fusioned,breaks = seq(1,LEN_GEN+BIN,by = BIN) , 
          xlab = "Position along the main genome (in bp)",
          ylab = "Number of contacts with plasmid pGM2",
          main =  paste("Contacts between plasmid and main genome with bin of",BIN,"bp") ,
          col = "deepskyblue4", 
          #      xlim =c(1e+06,2e+06)
)

int1 = h1$counts
int2 = h2$counts
#---------------------------------------------------------------------------------------------------
v=vector()
for(i in c(1:(length(h1$breaks)-1)) ) {v[i]=(h1$breaks[i]+h1$breaks[i+1])/2;}
#d=data.frame(v,h$counts);

# To calculate ratios: 
signal1=int1/cover1
signal2=int2/cover2


#--------------------------------------------------------------
# Plots:
#--------------------------------------------------------------
#  Plot of ratio of two signals:
plot(v,signal2,type="l",
     col="red",
     lwd=2.,
     xlab="Position along the Escherichia Coli genome",
     ylab="Interactions with plasmid / General cover ",
     main = paste("Ratio Contacts between plasmid and General Cover with bin of",BIN,"bp"),
     #xlim= c(1500000,1982000),
     ylim= c(0,0.03 ),
     cex=0.8
)
points(v,signal1,col="deepskyblue4",type="l",lwd=3)

#mats = c(1057964,1203180,1256343,1326562,1376485,1406840,1418870,1435317,1464772,1508258,
#    1508298,1530824,1560204,1636840,1684808)

mats = read.table("/home/axel/Bureau/data/vick/exp1/matS_cell.dat3")
points(mats$V1,mats$V1/mats$V1*.0,pch=19,col ="red")

legend("topright", inset=.082, 
       c("plasmid PGM2","plasmid PGM2-5xMatS"), 
       lwd=2.,
       lty=c(1,1),
       col=c("deepskyblue4","red"),
       bty = "n")

legend("topright", inset=.18, 
       c("MatS Sites"), 
       pch=19,
       col=c("red"),
       bty = "n")


#---------------
#  PLot of average of aligned signal around MatS: 
signal1[is.nan(signal1)] = 0
signal2[is.nan(signal2)] = 0
signal1[is.na(signal1)] = 0
signal2[is.na(signal2)] = 0

nbins =10   #  number of bins to look around the mats site
v1= matrix(0,nrow=length(mats$V1), ncol=nbins*2 +1)
v2= matrix(0,nrow=length(mats$V1), ncol=nbins*2 +1)
nb1= rep(0,nbins*2 +1)
nb2= rep(0,nbins*2 +1)

m=0
for(s in mats$V1) 
{
  m=m+1  
  s= floor(s/BIN)
  print(c("->",s))
  ii=0
  for(b in (s-nbins):(s+nbins) )
  { 
    ii=ii+1
    if(signal1[b] > 0) 
    {v1[[m,ii]] = signal1[b]
     nb1[ii] = nb1[ii] + 1 }
    if(signal2[b] > 0) 
    {v2[[m,ii]] = signal2[b]
     nb2[ii] = nb2[ii] + 1 }
    print(signal2[b])
  }
}

#  PLots:
plot(seq(1,21), colMedians(v2) / nb1, type="l" , lwd = 6.0, col = "red",
     xlab="Distance to matS site",
     ylab="Contact with plasmid",
     ylim= c(0,0.001))

points(seq(1,21),colMedians(v1) / nb2,type="l", lwd = 6.0, col = "deepskyblue4",
        xlab="Distance to matS site",
        ylab="Contact with plasmid"
)

abline(v =11, untf = FALSE, col ="red", lwd = 5.)
grid(nx = 11)


