#! /bin/bash

# The input files are in ./Data230622 folder
echo "All sequences in input files" 

cd ~/swDev/mon/Data230622
cat MP-MTAP-3repeat-pooled.fastq | echo "MTAP (all): $((`wc -l` / 4))"
cat MP-HCT116-3repeat-pooled.fastq | echo "HCT116 (all): $((`wc -l` / 4))"

#Extract sequences with q_score > 30
echo "MTAP: extract sequences with q_score > 30"
python ../Scripts/extractq.py ./MP-MTAP-3repeat-pooled.fastq ./MTAP_clusters.txt

echo "HCT116: extract sequences with q_score > 30"
python ../Scripts/extractq.py ./MP-HCT116-3repeat-pooled.fastq ./HCT116_clusters.txt

echo "Cut 18 chars from MTAP.txt"
cut -c 53-91 MTAP_clusters.txt | sed 's/./&:/12' | grep 'GTTGCTAGCAAT:'| sed 's/.*://' > MTAP_cut18.txt

echo "Cut 18 chars from HCT116.txt"
cut -c 53-91 HCT116_clusters.txt | sed 's/./&:/12' | grep 'GTTGCTAGCAAT:'| sed 's/.*://' > HCT116_cut18.txt

echo "Extract MTAP peptides"
python ../Scripts/translate.py ./MTAP_cut18.txt| grep -E 'CLS$'|grep -E '^S' |grep -v '*' |sed 's/.\{3\}$//' > MTAP_peptide.txt

echo "Extract HCT116 peptides"
python ../Scripts/translate.py ./HCT116_cut18.txt| grep -E 'CLS$'|grep -E '^S' |grep -v '*' |sed 's/.\{3\}$//' > HCT116_peptide.txt

echo "Find intersection between MTAP and HCT116 peptides"
python ../Scripts/set_v3.py MTAP_peptide.txt HCT116_peptide.txt
