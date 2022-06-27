# -*- coding: utf-8 -*-
"""

pip install XlsxWriter
pip install pandas
pip install bio
"""
#!/usr/bin/env python
# coding: utf-8

# ### Import used packages

from tkinter import TRUE
import xlsxwriter
import numpy as np
import matplotlib.pyplot as plt 
import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

# Main function 
def main():
    
    input_fastq = sys.argv[1]
    output_fastq = sys.argv[2]
    
    # Check if output file exists
    os.makedirs(os.path.dirname(output_fastq), exist_ok=True)
     
    # Strip out sequences with q_score>=30
    cluster = []
    for entry in SeqIO.parse(input_fastq, "fastq"): 
        if (sum(entry.letter_annotations['phred_quality'])/len(entry) >= 30):
            cluster.append(str(entry.seq))
    print("Clusters in: ", output_fastq, ':', len(cluster))
    
    with open(output_fastq, "w") as f:
        for sequence in cluster:
            f.write('{}\n'.format(sequence))

            
# python extractq.py inputFile outputFile

# Main program call
if __name__ == "__main__":
    main()
