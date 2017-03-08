# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 15:37:04 2014

@author: A.C
"""

# To calculate the GC content 

from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC_skew
import numpy as np


#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/bacillus_BS168/chr1.fa"), "fasta");
#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/bacillus_PY79/chr1.fa"), "fasta");
#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/bacillus_PY79_SMC_mutant/chr1.fa"), "fasta");
#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/bacillus_PY79_spoJ_mutant/chr1.fa"), "fasta");
#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/caulobacter/chr1.fa"), "fasta");
#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/vibrio_mazel/both_chrs.fa"), "fasta");
#record = SeqIO.read(open("/home/axel/Bureau/python/fasta/sacCer3_synth3/chr2.fa"), "fasta") ;  

def gc_signal(bin):
    record = SeqIO.read(open("/home/axel/Bureau/python/fasta/ecoli/chr1.fa"), "fasta");
    l = len(record.seq);
    n1 = int(l / bin)+1;
    
    gc_signal = np.zeros((n1, 1));
    gc_skew = np.zeros((n1, 1));

    for i in range(0,n1) :
        seq = record.seq[i*bin : i*bin + bin ];
        gc_signal[i]  =  GC(seq);
        gc_skew[i]  =  GC_skew(seq,window=bin);
        
    return gc_signal

#pos_gc = np.concatenate((pos,gc_signal ), axis=1);