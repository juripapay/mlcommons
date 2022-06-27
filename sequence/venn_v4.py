# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:47:08 2021

@author: lmw1e19
"""

# installation
# pip install matplotlib-venn

from collections import Counter
import matplotlib.pyplot as plt 
import os
from matplotlib_venn import venn3, venn3_circles,venn2,venn2_circles


# title
title = 'MTAP and HCT116'

#names = ['MTAP -i','MTAP +i']
names = ['MTAP','HCT116']
name1= names[0]
name2= names[1]

# folder
#folder = '/home/pearl061/swDev/mon/AnalysisTest3105'
folder = '/home/pearl061/swDev/mon/Analysis2206_second'

def text(sequences,name,out):
    """
    Function text creates a simple textfile where one sequence is written below the other.
    """   
    writer=open('%s-%s.txt'%(out,name),'w')
    for i in sequences:
        writer.write('{}\n'.format(i))
    writer.close()

def compare_datasets2(title,folder,apath,a,bpath,b,x=10):
    """
    This function compares the two datasets a and b. put in the path and name
    It also takes the simulated library into account.
    It generates text files that contain the respective overlapping or missing sequences.
    """
    
    outfolder=folder+'/compare%s/'%(title)
    out=outfolder+'cCX5'
    if not os.path.exists(outfolder):
             os.makedirs(outfolder)
    f=open(apath,'r')
    pcr=f.readlines()
    pcr=[seq[:6] for seq in pcr]
    pcrCount=Counter(pcr)
    pcr=[i for i in pcrCount if pcrCount[i]>x]
    f.close()

    f=open(bpath,'r')
    pre=f.readlines()
    pre=[seq[:6] for seq in pre]
    preCount=Counter(pre)
    pre = [i for i in preCount if preCount[i]>x]
    f.close()

    #make all the lists
    pcrpre=pcr+pre    
    pcrpreCount=Counter(pcrpre)
    
    #build the intersections based on sequences
    onlypcr=[i for i in pcr if preCount[i]<=x]
    onlypre=[i for i in pre if pcrCount[i]<=x]
    inall=[i for i in list(pcrpreCount) if pcrCount[i]>x and preCount[i]>x]

    #generate lists with the sequences of each group 
    text(onlypcr,'%sonly%s'%(x,a),out)
    text(onlypre,'%sonly%s'%(x,b),out)
    text(inall,'%sin all'%(x),out)   

    sets=Counter()
    sets['10']=len(onlypcr)
    sets['01']=len(onlypre)
    sets['11']=len(inall)
    
    fig=plt.figure()
    v=venn2(subsets = sets, set_labels = ('%s'%(a),'%s'%(b)))
    c=venn2_circles(subsets = sets, linestyle='-', linewidth=1, color="black")

    h,l=[],[]
    for i in sets:
        h.append(v.get_patch_by_id(i))
        l.append(sets[i])
        if sets[i]>0:
            v.get_label_by_id(i).set_text("")
    # plt.title('Drop out screen 1 and 2')
    plt.legend(title='unique sequences',handles=h, labels=l,bbox_to_anchor=(0.8,0.2), loc="upper left")
    plt.text(1,0.7,'total A:%s, total B:%s'%(sum(pcrCount.values()), sum(preCount.values())))
    plt.title(title+' over %s'%(x))
    plt.savefig("%s-venn-5.png" %(out))
    plt.show()


def compare_datasets(title,folder,apath,a,bpath,b):
    """
    This function compares the two datasets a and b. put in the path and name
    It also takes the simulated library into account.
    It generates text files that contain the respective overlapping or missing sequences.
    """
    
    outfolder=folder+'/compare%s/'%(title)
    out=outfolder+'cCX5'
    if not os.path.exists(outfolder):
             os.makedirs(outfolder)
    f=open(apath,'r')
    pcr=f.readlines()
    pcr=[seq[:6] for seq in pcr]
    pcrCount=Counter(pcr)
    pcr=[i for i in pcrCount]
    f.close()
    f=open(bpath,'r')
    pre=f.readlines()
    pre=[seq[:6] for seq in pre]
    preCount=Counter(pre)
    pre = [i for i in preCount]
    f.close()
    #make all the lists
    pcrpre=pcr+pre    
    pcrpreCount=Counter(pcrpre)
    
    #build the intersections based on sequences
    onlypcr=[i for i in pcr if preCount[i]==0]
    onlypre=[i for i in pre if pcrCount[i]==0]
    inall=[i for i in list(pcrpreCount) if pcrCount[i]>0 and preCount[i]>0]
    #generate lists with the sequences of each group 
    text(onlypcr,'only%s'%(a),out)
    text(onlypre,'only%s'%(b),out)
    text(inall,'in all',out)    
    sets=Counter()
    sets['10']=len(onlypcr)
    sets['01']=len(onlypre)
    sets['11']=len(inall)
    fig=plt.figure()
    v=venn2(subsets = sets, set_labels = ('%s'%(a),'%s'%(b)))
    c=venn2_circles(subsets = sets, linestyle='-', linewidth=1, color="black")

    h,l=[],[]
    for i in sets:
        h.append(v.get_patch_by_id(i))
        l.append(sets[i])
        if sets[i]>0:
            v.get_label_by_id(i).set_text("")
    # plt.title('Drop out screen 1 and 2')
    #plt.legend(title='unique sequences',handles=h, labels=l,bbox_to_anchor=(1,0.27), loc="upper left")

    plt.legend(title='unique sequences',handles=h, labels=l,bbox_to_anchor=(0.8,0.2), loc="upper left")
    
    plt.text(1,0.7,'total A:%s, total B:%s'%(sum(pcrCount.values()), sum(preCount.values())))
    plt.title(title)
    plt.savefig("%s-venn.png" %(out))
    plt.show()


# Main program
do1=folder+'/'+name1+'/'+name1+'_cSX5.txt'
do2=folder+'/'+name2+'/'+name2+'_cSX5.txt'
print("do1: ", do1)
print("do2: ", do2)

compare_datasets(title,folder, do1,name1 , do2, name2)   
compare_datasets2(title,folder, do1,name1 , do2, name2,x=5)   
compare_datasets2(title,folder, do1,name1 , do2, name2,x=10)  
