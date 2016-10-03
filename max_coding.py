# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:24:02 2016

@author: daxinggao
"""
fh=open("/Users/daxinggao/dna2.fasta")
fh.seek(0)
dic={}
for line in fh:
    if line.startswith(">"):
        name=line.rstrip()
        dic[name]=''
    else:
        dic[name]=dic[name]+line.rstrip()
lenth={}
for k,v in dic.iteritems():
    lenth[k]=len(v)
mini=sorted(lenth.values())[0]
maxi=sorted(lenth.values())[-1]
for key in lenth:
    if lenth[key]==mini:
        print "min length",key
    if lenth[key]==maxi:
        print "max length",key

def rffinder(dna,n):
    pos=-1
    maxx=0
    for i in range(n,len(dna),3):
        if dna[i:i+3]=="ATG":
            pos=i
            for j in range(i+3,len(dna),3):
                if dna[j:j+3]=="TAA" or dna[j:j+3]=="TAG" or dna[j:j+3]=="TGA":
                    if j+3-pos>maxx:
                        maxx=j+3-pos
                        posmax=pos
                    break
    #print maxx
    #print posmax
    return maxx
    
def findrep(n):
    maxcount=0
    lis=[]
    for k,v in dic.iteritems():
        for i in range(len(v)-n+1):
            count=0
            for k2,v2 in dic.iteritems():
                count=count+v2.count(v[i:i+n])
            if count>maxcount:
                maxcount=count
                maxseq=v[i:i+n]
            if count==5:
                lis.append(v[i:i+n])
            
    print maxcount
    print maxseq
    return lis
    
    
for k,v in dic.iteritems():
    print rffinder(v,1)

for k in dic.keys():
    if re.search('gi\|142022655\|gb\|EQ086233\.1\|16',k):
        print k
        
for k,v in dic.iteritems():
    print findrep(v,6)    
            
    
    
    