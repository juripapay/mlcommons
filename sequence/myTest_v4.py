# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:29:00 2022

@author: lmw1e19
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
import os, yaml, sys, argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

# ## Text writing functions
# These functions are later used for analysis by the program itself as well as an input for the database
# ### Plain text

# This function writes a plain text file of a list of sequences.
def plain_text(sequences,name,output_folder = ''):
    writer = open('%s%s.txt'%(output_folder, name),'w')
    for sequence in sequences:
        writer.write('{}\n'.format(sequence))
    writer.close()
    print('plain text for ', name, ' saved as %s%s.txt'%(output_folder, name))

# ### Tabulated text from Counter
#This function writes a text file from a counter, which is easily readable for the user.
def counter_text(sequence_count ,name, output_folder = ''):
    # Generate a list that sorts the counter from most to least common.
    most_common_seqs = sequence_count.most_common()
    # Defines the number of total sequences
    total_seqs = sum(sequence_count.values())
    writer = open('%s%s-mc.txt'%(output_folder, name), 'w')
    for sequence, occurence in most_common_seqs:
        # Write the sequence of the peptide, the occurence and the fraction of this peptide in the sequenceing run.
        writer.write('{}  {}  {:1.4f}\n'.format(sequence, occurence, occurence/total_seqs*100))
    writer.close()
    print('counter text file for ', name, ' saved as %s%s-mc.txt'%(output_folder, name))

# ### The run analyzer main function
# This function has a lot of tasks
def run_analyzer(primer_dict, clusters, name, output_folder = '',
                forward_primer = 'Cfa-000', starting_aminoacid = 'S'):
    # Define a new output folder path for all the primary data analysis.
    output_primary = output_folder
    # Generate than new path if it doesn't exist already.
    if not os.path.exists(output_primary):
        os.makedirs(output_primary)
    # Check which if a correct fw primer was used to assign the starting point for data processing 
    # and NheI restriction site surrounding sequence as well as the N-intein sequence and cropping number.
    if forward_primer in primer_dict:
        start = primer_dict[forward_primer][0]
        nheI_restriction = primer_dict[forward_primer][1]
        test_length = len(nheI_restriction)
        n_intein = primer_dict[forward_primer][3]
        cropping = primer_dict[forward_primer][4]
    else:
        # Exit the program so the right primer can be used
        print('please insert a valid primer as Int-123')
        return
    # Crop the sequences if they contain the right NheI secuence at the right position to start after the 
    # NheI restriction site and remove any overhang that would disrupt later translation
    # Copy the sequences when you find the nheI_restriction EXACTLY at 52nd position. We can collect stats on other positions.

    intein_dna = [seq[start+test_length:-cropping] for seq in clusters if 
                  seq[start:start+test_length] == nheI_restriction]
    
    restriction_positions = []
    for item in clusters:
        itemStr = str(item)
        restriction_positions.append(itemStr.find(nheI_restriction))
    
    print("Total number of clusters in: ", name, len(clusters))
    print("Restriction enzyme found in the right position: ", name, len(intein_dna))
    print("Restriction enzyme found in ANY position: ", name, Counter(restriction_positions))

    # Create a list for the extein DNAs
    extein_dna = []
    # Create a list for the extein amino acids (peptides) starting with C (starting amino acid) and other.
    extein_S = []
    # Fill the lists by iterating through the intein DNA
    for sequence in intein_dna:
        # Convert the string sequence into a biopython Seq object and translate it to amino acid sequence
        # and transform it to a string again for easier data storage
        seq = str(Seq(sequence).translate())
    
        # Check if the translated sequence is in line with the expected sequence for the N-intein ('CLS'), the sequence must also contain CLS at the right position.
        if seq[6:6+len(n_intein)]== n_intein:
            # Add the peptides with the desired amino acid to the respective cysteine or other lists for peptides.
            # The extein_dna list will contain only the exteins of the desired starting amino acid.
            if seq.startswith(starting_aminoacid):
                extein_S.append(seq[:6])
                extein_dna.append(sequence[:18])

    # Remove all stopcodons to obtain only the cyclic ones and save as a plain text file.
    print (name, ",Number of exteins:", len(extein_S))
    # Cyclic peptidea are erxteins without *, X
    cyclic_peptides = [seq for seq in extein_S if '*' not in seq and 'X' not in seq]
    plain_text(cyclic_peptides, name+'_c'+starting_aminoacid+'X5', output_folder = output_primary)
    print(name, ", Number of cyclic peptides: ", len(cyclic_peptides))
      
    # intein contains restriction about 80 chars, do we need this file?
    intein_dna_count = Counter(intein_dna)
    counter_text(intein_dna_count, name+'_intein', output_folder = output_primary) 
    
    # Count the occurence of exteins (may contain *,X), do we need it?
    extein_dna_count = Counter(extein_dna)
    counter_text(extein_dna_count, name+'_exteinDNA', output_folder = output_primary)
    
    # Count frequence of exteins AA, do we need?
    extein_S_count = Counter(extein_S)
    counter_text(extein_S_count, name+'_'+starting_aminoacid+'-extein', output_folder = output_primary)

    # Frequency of cyclic peptides
    cyclic_peptides_count = Counter(cyclic_peptides)
    counter_text(cyclic_peptides_count, name+'_cSX5',output_folder = output_primary)
    
    # Generate an excel workbook to store the absolute and unique copy numbers.
    # How to interpret the sprwadsheet?
    wb = xlsxwriter.Workbook("%s%s_stats.xlsx"%(output_primary,name))
    # Add a worksheet to the workbook.
    ws = wb.add_worksheet('summary')
    # Write tha data in the worksheet.
    ws.write(0, 0, name)
    ws.write(1, 0, 'total clusters') 
    ws.write(2, 0, 'total exteins') 
    ws.write(3, 0, 'unique exteins')
    ws.write(4, 0, 'total exteins pep')
    ws.write(5, 0, 'unique exteins pep')
    ws.write(6, 0, 'total CP')
    ws.write(7, 0, 'unique CP')
    ws.write(9, 0, 'other total')
    ws.write(10, 0, 'other unique')
    ws.write(1, 1, len(clusters)) # all sequences, with a condition q_score>20
    ws.write(2, 1, len(extein_dna)) # sequence contains restriction and DNA equivalent of CLS, 18 chars
    ws.write(3, 1, len(extein_dna_count)) # frequency
    ws.write(4, 1, len(extein_S)) # AA sequence
    ws.write(5, 1, len(extein_S_count)) # frequency AA seqquence
    ws.write(6, 1, len(cyclic_peptides)) # peptides which don't contain *, X
    ws.write(7, 1, len(cyclic_peptides_count)) # frequency of peptides  
    wb.close()
    return(cyclic_peptides, extein_dna)

# # Numeric analysis of the runs
# ## DNA plot
# The baseplot generates a plot to show the base distribution per position
# Use a 2D array, it will be simpler
def baseplot(extein_dna, name, output_folder = ''):
    '''
    The baseplot function calculates the percentage base per position and plots it.
    '''
    # Define the DNA bases
    bases=['T','C','A','G']
    # The length of the extein is 18 bp
    length=18
    # Generate a numpy array filled with zeros for every base and every position
    basecount=np.zeros((length,4))
    # Iterate through all the exteins.
    for sequence in extein_dna:
        # Iterate through the positions and bases of the extein.
        for position, base in enumerate(sequence):
            # Identify the index of the base in the bases list and add +1 in the respective position of the matrix.
            for i in [i for i,x in enumerate(bases) if x==base]:
                basecount[position][i]+=1
    # Normalize the matrix to get 100 % for each position
    for position in range(length):
        basecount[position,:]*=100/sum(basecount[position,:])
    # Set the x-axis numbers
    x=np.arange(1,length+1)
    # Set the colors for the bases
    colors=['red','blue','green','black']
    # Add grey lines at 25, 50 and 75 %
    plt.axhline(y=50,color='dimgray',linewidth=0.5)
    plt.axhline(y=25,color='dimgray',linewidth=0.5)
    plt.axhline(y=75,color='dimgray',linewidth=0.5)
    # Draw the dots for all four bases
    for i in range(4):
         plt.scatter(x=x, y=basecount[:,i], label=bases[i] , color= colors[i],  
                     marker='.', s=30)
    # x-axis label
    plt.xlabel('Position') 
    plt.xticks(np.arange(1,length+1,1))
    # frequency label 
    plt.ylabel('% Base') 
    plt.yticks(np.arange(0,101,10))
    # plot title 
    plt.title(name)
    # showing legend 
    plt.legend(loc='upper right')
    # save the plot as pdf
    plt.savefig("%s%s-bases.pdf" %(output_folder,name))
    # function to show the plot 
    plt.show()

# ## Plot the copy number
# ### Calculate the plot function
def counter_count(sequences, output_folder = '', name = 'xy', output = 'yes'):
    # Generate a counter from the sequences
    seq_counter = Counter(sequences)
    
    # The number of times a copy number is present is counted here.
    copy_numbers = []
    # Iterate through the counter object by ranking it.
    for entries in seq_counter.most_common():
        # Add the copy numbers (second element of the 'entries')
        copy_numbers.append(entries[1])
    # Count how many times each copy number occurs
    copy_numbers_counter = Counter(copy_numbers)
    
    # Transform into a list of tupels:
    copy_numbers_xy = copy_numbers_counter.most_common()
    # Extract into two lists to generate the x and y values for the plot.
    x = [i[0] for i in copy_numbers_xy]
    y = [i[1] for i in copy_numbers_xy]
    if output == 'yes':
        dataframe = pd.DataFrame(zip(x,y), columns = ['copy number', 'occ of copy number'])
        dataframe.to_csv('%s%s.csv'%(output_folder, name))
    return (x, y)

# ### Draw the copy number plot
def counter_count_plot(x_list, y_list, namelist, name, output_folder = ''):
    # Define the colours for the different datasets (max = 4)
    colors = ['b', 'r', 'darkorange', 'g']
    # Plot up to 4 datasets
    for index, x in enumerate(x_list):
        # Check that there are not too many data sets
        if index > 3:
            print ('Only four data sets allowed')
            return
        else:
            # plot the datasets individually
            y = y_list[index]
            label = namelist[index]
            plt.scatter(x = x, y = y, label = label, color = colors[index], marker = 'o', s = 5)
    # Add the labels
    plt.xlabel('Copy number (# of times unique sequence was found)')
    plt.ylabel('Frequency')
    plt.title(name)
    # Scale the axes to log
    plt.xscale('log')
    plt.yscale('log')
    # Save the plot
    plt.savefig('%s%soccurence.pdf'%(output_folder, name))
    plt.show()

#
# Main function
# python myTest_v4.py --config ./sequenceConfig_2306.yaml
#
def main():
 
     # Read command line arguments
    parser = argparse.ArgumentParser(description='CloudMask command line arguments',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--config', default=os.path.expanduser('./sequenceConfig.yaml'), help='path to config file')
    command_line_args = parser.parse_args()

    configFile = os.path.expanduser(command_line_args.config)
    
     # Read config file  
    with open(configFile, 'r') as stream:
        args = yaml.safe_load(stream)
    mtapFile =  os.path.expanduser(args['mtapdata'])
    hct116File = os.path.expanduser(args['hct116data'])
    input_fastq = []
    input_fastq.append(mtapFile)
    input_fastq.append(hct116File)
    outfolder = os.path.expanduser(args['outfolder'])
    
    mtapName = args['mtapName']
    hct116Name = args['hct116Name']
    fastq_names = [mtapName, hct116Name]
    print("input_fastq: ", input_fastq)
    print("fastq_names: ", fastq_names)
   
    # Define the folder where the data is meant to be stored by the end of the analysis. 
    # This folder doesn't need to exist yet. All analysis will go there.
    # Some sub-folders will also be created
    # Create MTAP and HCT116 folders 
    output_folders = [outfolder + i + '/' for i in fastq_names] 
    #
    # Set the primer that was used in this run
    forward_primer= 'Cfa-000'

    # Set the starting amino acid (cysteine C is default)
    starting_aminoacid = 'S'

        # ## Main function for primary fastq processing
    # The run analyzer should be run for every new run. The individual functions can be called at a later point as well.
    
    # ### Define the allowed primers
    # For every intein there is at least 1 primer, which determines the starting point of sequencing. There are some inteins, where more than one is present. Here, the primers are defined and the starting point is given.

    # Generates the primer dictionary
    # Every primer has a name (3 letters for the intein - 3 numbers for the ID, 
    # a distance (in bp) from the start of the sequencing to the start of the NheI site,
    # the intein-specific sequence for the splice junction (NheI surrounding bases),
    # the sequence of the primer itself (this is for checking if the right primer was used),
    # the N-intein sequence that is expected after the extein (in amino acid),
    # and a correction number to crop the sequence so that it results in triplets for translation 
    # (the overhang bases when having a 150 bp sequencing read).
    primer_dict = {'Cfa-000': [52, 'GTTGCTAGCAAT',
                            'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTATGACATAGGCGTCGGTGAGCCAC',
                            'CLS', 1]}
    # Check if all the primers have been assigned
    for keys, items in primer_dict.items():
        print(keys)
 

    # ## Load the clusters from the run and call the main analyzer function and assign the outputs
    # Here, all clusters with a quality over 30 on average are allowed.
    # The parameters for this were already set on the top of the page.

    # Create empty lists to store the results of the main function
    # These lists are used by future functions and store the cyclic peptide sequences and their coding DNA
    cyclic_peptides_byrun = []
    exteins_byrun = []

    # This operation is built to analyze a few runs in a row
    for path, name, output_folder in zip(input_fastq, fastq_names, output_folders):
        # The individual clusters are extracted as strings from the fastq file that was given at the start of the sheet.
        # Only entries with q-scores over 30 are allowed.

        clusters=[str(entry.seq) for entry in SeqIO.parse(path, "fastq") 
            if (sum(entry.letter_annotations['phred_quality']))/len(entry)>=30]
        
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # How many times the sequence occurs. How many unique sequences?
        clusterc=Counter(clusters)
        # How many unique sequences? 
        # len(clusters)=4119189, len(clusterc)=1435689

        # Save in a file
        counter_text(clusterc, name+'_intein_all', output_folder = output_folder) 
        
        # mc for most common first 50 characters
        # Why do we need it of these 50 characters are cut off anyway
        '''
        clusters2=[i[:50] for i in clusters]
        clusterc2=Counter(clusters2)
        counter_text(clusterc2, name+'_intein_crop_all', output_folder = output_folder) # mc for most common
        '''

        # Print a control if the correct dataset was loaded and how many clusters there are.
        # This gives an idea how long the analysis might take
        print(name, 'loaded', len(clusters), ' total clusters')
        print(name, 'loaded', len(clusterc), ' number of unique clusters (150 characters)')
        
        # The main run analyzer function is called for that specific run
        cCX5, extein_dna = run_analyzer(primer_dict, clusters, name, output_folder = output_folder, 
                                        forward_primer = forward_primer,
                                        starting_aminoacid = starting_aminoacid)
        
        # Peptides array, Peptide:SVLSLF                      
        cyclic_peptides_byrun.append(cCX5)
        # Corresponding DNA array, DNA:AGCGTCTTGTCCCTCTTT     
        exteins_byrun.append(extein_dna)
        print(name, ' saved')

        # ## Execute the numeric analysis

    # ### The DNA plot
    # For all data sets that were generated the same name and folder structure is used
   
    for index, extein_dna in enumerate(exteins_byrun):
        name = fastq_names[index]
        output_folder = output_folders[index]
        baseplot(extein_dna, name, output_folder=output_folder)
    
    # ### The copy number plot
    # For all data sets that were generated the same name and folder structure is used.
    # Individual runs can be plotted as well as up to 4 different runs. For this, the x and y are saved as CSV.
    x_list = []
    y_list = []
    name_list = []
    
    for index, cyclic_peptides in enumerate(cyclic_peptides_byrun):
        name = fastq_names[index]
        output_folder = output_folders[index]
        # Calculate the copy numbers and their occurences
        x, y = counter_count(cyclic_peptides, output_folder = output_folder, name = name)
        # Plot the data as lists
        counter_count_plot([x], [y], [name],'', output_folder = output_folder)
        # Add the data to the list for easy plotting of multiple datasets
        x_list.append(x)
        y_list.append(y)
        name_list.append(name)

# python OpenSX5_distribution3105.py

# Main program call
if __name__ == "__main__":
    main()
