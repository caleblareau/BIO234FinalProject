# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
import random 

fl = 'sim_600_t0'
hapFile = '../input/'+fl+'.hap'
phenoFile = '../input/'+fl+'.sample'

#Read in allele data
raw = []
with open(hapFile,'r') as f:
   for line in f:
    #only read in about 1% of the variants
    if(random.uniform(0,1) < 0.01):
      raw.append([ int (x) for x in line.split('\t') ])

#Read in ethnicity data
pheno = []
with open(phenoFile,'r') as f: 
    for line in f:
       pheno.append(line.split()[0])

G = np.matrix(raw)

#Number of variants per person per position (0,1,2)
geneticMatrix = np.zeros((G.shape[0],G.shape[1]/2))
for i in range(0,G.shape[1],2):
    geneticMatrix[:,] = (G[:,i]+G[:,i+1])
    
#Covariance Matrix
genCov = np.cov(geneticMatrix, rowvar = 0)

#Run PCA
results = PCA(genCov)

#assign colors for PCA plot
colors_dict = {}
pheno_keys = list(set(pheno))
for i in range(len(pheno_keys)):
    colors_dict[pheno_keys[i]] = i

colors = []
for i in range(len(pheno)):
    colors.append(colors_dict[pheno[i]])

#plot first two components of the PCA
plt.scatter(results.Y[:,0],results.Y[:,1], c = colors)
plt.ylabel('PC1')
plt.xlabel('PC2')
pltfile = '../output/'+fl+'.png'
plt.savefig(pltfile)

