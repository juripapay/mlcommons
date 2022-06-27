# Translate to AA sequence
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Read file into array
fileName = sys.argv[1]
my_file = open(fileName, "r")
  
# reading the file
data = my_file.read()
  
# replacing end splitting the text 
# when newline ('\n') is seen.
data_into_list = data.split("\n")
my_file.close()

seqAA = []
for item in data_into_list:
	seq = str(Seq(item).translate())
	seqAA.append(seq)
	print(seq)

