---
jupyter:
  jupytext:
    formats: ipynb,md,py
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

### Week 2 Project
My idea for this week's project contribution is to take 10 covid-19 genome sequences from the beginning of the outbreak in the U.S. 
and 10 more recent covid-19 sequences and compare them against a reference spike protein to see if the virus has mutated 
at all since the time it's been spreading here in the country.


### Combining Files
All the genomes are downloaded into separate fasta files, so the first thing to do is combine all the genomes into separate files
depending on if their recent sequences or early sequences.

```python
import os

recentPath = "./covidSeqs/recentSeqs/"

with open(recentPath + "recentSequences.fasta", 'w') as newfile:
    for genomeFile in os.listdir(recentPath):
        if genomeFile != "recentSequences.fasta":
            file = open(recentPath + genomeFile, 'r')
            genome = file.read()
            newfile.write(genome + "\n")
```

```python
### we do the same thing for the early sequences

earlyPath = "./covidSeqs/earlySeqs/"

with open(earlyPath + "earlySequences.fasta", 'w') as newfile:
    for genomeFile in os.listdir(earlyPath):
        if genomeFile != "earlySequences.fasta":   
            file = open(earlyPath + genomeFile, 'r')
            genome = file.read()
            newfile.write(genome + "\n")
```

### Virus Alignment
I'm using virulign for aligning the viruses against a reference protein. The repository for the virulign program contains the data
for the reference protein. Below, I have downloaded both the program and the repository.

```python
!/Users/ericinman/bin/virulign
```

I also downloaded the tutorials and the program repository.

```python
!git clone https://github.com/rega-cev/virulign ../../../virulign
```

### Alignment
Now we need to align the genomes with the reference protein so that we can generate a distance matrix.

We're going to do this for both the recent sequences and the early sequences.

```python
!/Users/ericinman/bin/virulign ../../../virulign/references/SARS-CoV-2/S.xml ./covidSeqs/recentSeqs/recentSequences.fasta \
--exportAlphabet Nucleotides --exportKind PositionTable > ./covidSeqs/recentSeqs/recent_position_table.csv
```

```python
!/Users/ericinman/bin/virulign ../../../virulign/references/SARS-CoV-2/S.xml ./covidSeqs/earlySeqs/earlySequences.fasta \
--exportAlphabet Nucleotides --exportKind PositionTable > ./covidSeqs/earlySeqs/early_position_table.csv
```

### Reading the data into a pandas dataframe

```python
import pandas as pd
recent_position_table = pd.read_csv('./covidSeqs/recentSeqs/recent_position_table.csv')
early_position_table = pd.read_csv('./covidSeqs/earlySeqs/early_position_table.csv')
```

These are our tables below:

```python
recent_position_table
```

```python
early_position_table
```

### Combining Charts and Taking Consensus Sequences

```python
# the early sequences will be on top, recent sequences on the bottom
combinedTables = pd.concat([early_position_table, recent_position_table], axis = 0)
```

Now we take the consensus sequences from all three of the tables we've created.

```python
early_concensus_seq = early_position_table.drop('seqid',axis=1).mode(axis=0).T[0]
early_concensus_seq
```

```python
recent_concensus_seq = recent_position_table.drop('seqid', axis = 1).mode(axis = 0).T[0]
recent_concensus_seq
```

```python
combined_concensus_seq = combinedTables.drop('seqid', axis = 1).mode(axis = 0).T[0]
combined_concensus_seq
```

```python
early_position_table = early_position_table.set_index('seqid')
recent_position_table = recent_position_table.set_index('seqid')
combinedTables = combinedTables.set_index('seqid')
```

### Determine which samples are farthest from the concensus sequence

```python
early_distance_from_concensus_seq = early_position_table.apply(lambda row: sum(row != early_concensus_seq),axis=1)
early_distance_from_concensus_seq_sorted = early_distance_from_concensus_seq.sort_values(ascending=False)
early_distance_from_concensus_seq_sorted
```

```python
recent_distance_from_concensus_seq = recent_position_table.apply(lambda row: sum(row != recent_concensus_seq),axis=1)
recent_distance_from_concensus_seq_sorted = recent_distance_from_concensus_seq.sort_values(ascending=False)
recent_distance_from_concensus_seq_sorted
```

```python
combined_distance_from_concensus_seq = combinedTables.apply(lambda row: sum(row != combined_concensus_seq),axis=1)
combined_distance_from_concensus_seq_sorted = combined_distance_from_concensus_seq.sort_values(ascending=False)
combined_distance_from_concensus_seq_sorted
```

### Constructing a distance matrices for our sequences

```python
earlySequids = early_distance_from_concensus_seq_sorted.index
recentSequids = recent_distance_from_concensus_seq_sorted.index
combinedSequids = combined_distance_from_concensus_seq_sorted.index
earlyDistances = {}

for i,seqid1 in enumerate(earlySequids):
    earlyDistances[seqid1,seqid1]=0
    for j in range(i+1,len(earlySequids)):
        seqid2 = earlySequids[j]
        earlyDistances[seqid1,seqid2] = sum(early_position_table.loc[seqid1] != early_position_table.loc[seqid2])
        earlyDistances[seqid2,seqid1] = earlyDistances[seqid1,seqid2]
earlyDistances = pd.Series(earlyDistances).unstack()
earlyDistances
```

```python
recentDistances = {}

for i,seqid1 in enumerate(recentSequids):
    recentDistances[seqid1,seqid1]=0
    for j in range(i+1,len(recentSequids)):
        seqid2 = recentSequids[j]
        recentDistances[seqid1,seqid2] = sum(recent_position_table.loc[seqid1] != recent_position_table.loc[seqid2])
        recentDistances[seqid2,seqid1] = recentDistances[seqid1,seqid2]
recentDistances = pd.Series(recentDistances).unstack()
recentDistances
```

```python
combinedDistances = {}

for i,seqid1 in enumerate(combinedSequids):
    combinedDistances[seqid1,seqid1]=0
    for j in range(i+1,len(combinedSequids)):
        seqid2 = combinedSequids[j]
        combinedDistances[seqid1,seqid2] = sum(combinedTables.loc[seqid1] != combinedTables.loc[seqid2])
        combinedDistances[seqid2,seqid1] = combinedDistances[seqid1,seqid2]
combinedDistances = pd.Series(combinedDistances).unstack()
combinedDistances
```

### Utilizing biopython
Lastly, we'll use biopython to build our phylogenies

```python
from Bio.Phylo.TreeConstruction import DistanceMatrix
import numpy as np

earlyMatrix = np.tril(earlyDistances.values).tolist()
for i in range(len(earlyMatrix)):
    earlyMatrix[i] = earlyMatrix[i][:i+1]
earlydm = DistanceMatrix(list(earlyDistances.index), earlyMatrix)
```

```python
recentMatrix = np.tril(recentDistances.values).tolist()
for i in range(len(recentMatrix)):
    recentMatrix[i] = recentMatrix[i][:i+1]
recentdm = DistanceMatrix(list(recentDistances.index), recentMatrix)
```

```python
combinedMatrix = np.tril(combinedDistances.values).tolist()
for i in range(len(combinedMatrix)):
    combinedMatrix[i] = combinedMatrix[i][:i+1]
combineddm = DistanceMatrix(list(combinedDistances.index), combinedMatrix)
```

### Now construct our tree

```python
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()
earlyTree = constructor.nj(earlydm)
recentTree = constructor.nj(recentdm)
combinedTree = constructor.nj(combineddm)
```

### Now draw our trees

```python
%matplotlib inline

from Bio import Phylo
earlyTree.ladderize()   # Flip branches so deeper clades are displayed at top
recentTree.ladderize()
combinedTree.ladderize()
Phylo.draw(earlyTree)
Phylo.draw(recentTree)
Phylo.draw(combinedTree)
```

**Please see the guidance at the top of the page for what to try**

```python

```
