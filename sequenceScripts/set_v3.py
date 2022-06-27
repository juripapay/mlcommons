# Find the intersection between two peptide sets
# Run it from Data3105 folder
# python ../Scripts/set_v3.py MTAP_peptide.txt HCT116_peptide.txt
import sys
import matplotlib.pyplot as plt
from collections import Counter
# For 3 sets use venn3
from matplotlib_venn import venn2

def saveFile(fileName, setName):
    with open(fileName,'w') as f:	
        for item in setName:
            longStr = item + item + '\n'
            f.write(longStr)
            longStr = ''
    f.close()


# MTAP_peptide.txt
mtapFile = sys.argv[1]
# HCT116 peptides
hct116File = sys.argv[2]

# Read MTAP
mtapList = []
with open(mtapFile) as f:
    mtapList = f.read().splitlines()
mtapSet = set(mtapList)
countMTAP = Counter(mtapList)

# Make a copy of dictionary
countMTAP5 = countMTAP.copy()
countMTAP10 = countMTAP.copy()

# Filter out elements with frequency <5
for k, v in countMTAP5.copy().items():
    if v < 6:
        countMTAP5.pop(k)

# Filter out elements with frequency <10
for k, v in countMTAP10.copy().items():
    if v < 11:
        countMTAP10.pop(k)

# Make sets
mtapSet5 = set(countMTAP5.keys())
mtapSet10 = set(countMTAP10.keys())

# Read HCT116
hct116List = []
with open(hct116File) as f:
    hct116List = f.read().splitlines()
hct116Set = set(hct116List)
countHCT116 = Counter(hct116List)

# Make a copy of dictionary
countHCT116_5 = countHCT116.copy()
countHCT116_10 = countHCT116.copy()

# Filter out elements with frequency <5
for k, v in countHCT116_5.copy().items():
    if v < 6:
        countHCT116_5.pop(k)

# Filter out elements with frequency <10
for k, v in countHCT116_10.copy().items():
    if v < 11:
        countHCT116_10.pop(k)

# Make sets
hct116Set5 = set(countHCT116_5.keys())
hct116Set10 = set(countHCT116_10.keys())

# Draw venn diagram
venn2([mtapSet, hct116Set], ('MTAP', 'HCT116'))
plt.show()

# Peptides occuring >5 
venn2([mtapSet5, hct116Set5], ('MTAP5', 'HCT116_5'))
plt.show()

# Peptides occuring >10
venn2([mtapSet10, hct116Set10], ('MTAP10', 'HCT116_10'))
plt.show()

# We are interested in HCT116 dropouts
dropout_HCT116 = hct116Set.difference(mtapSet)
dropout_HCT116_5 = hct116Set5.difference(mtapSet5)
dropout_HCT116_10 = hct116Set10.difference(mtapSet10)

# Save dropouts in files
saveFile("dropout_HCT116.txt", dropout_HCT116)
saveFile("dropout_HCT116_5.txt", dropout_HCT116_5)
saveFile("dropout_HCT116_10.txt", dropout_HCT116_10)
