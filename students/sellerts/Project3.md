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

## Project 3


This week I continued my analysis from last week, looking further into the generated phylogenetic trees that I arrived at last time.



### Imports

```python
!pip install biopython
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import pandas as pd
```

### Loading the given data from the position table

```python
position_table = pd.read_csv('../../data/position_table.csv')
```

```python
results = position_table.describe()
results
```

```python
position_table
```

### Pull out the concensus sequence

```python
concensus_seq = position_table.drop('seqid',axis=1).mode(axis=0).T[0]
concensus_seq
```

```python
position_table = position_table.set_index('seqid')
```

### Determine which samples are farthest from the concensus sequence

```python
distance_from_concensus_seq = position_table.apply(lambda row: sum(row != concensus_seq),axis=1)
distance_from_concensus_seq_sorted = distance_from_concensus_seq.sort_values(ascending=False)
distance_from_concensus_seq_sorted
```

### Select Sets of sequences


Here I break from the outline, and I pull out four sets of distance data from the sorted list.
All the following analysis will be in order of decreasing distance from the concensus sequence, since
that is how the data is sorted.

I also decided to look at the data in groups of eight instead of ten, since that will allow us to see the graphs better

```python
subset_seqs = []
for i in range(0, 32, 8):
    subset_seqs.append(distance_from_concensus_seq_sorted[i:i + 8].index)
subset_seqs
```

### Calculate distances for each set of sequences

```python
def calc_distances(subset_seqs):
    distances = {}
    for i,seqid1 in enumerate(subset_seqs):
        distances[seqid1,seqid1]=0
        for j in range(i+1,len(subset_seqs)):
            seqid2 = subset_seqs[j]
            distances[seqid1,seqid2] = sum(position_table.loc[seqid1] != position_table.loc[seqid2])
            distances[seqid2,seqid1] = distances[seqid1,seqid2]
    distances = pd.Series(distances).unstack()
    return distances
    
distance_sets = [calc_distances(sequence) for sequence in subset_seqs]
distance_sets
```

As you can see, the distances from the concensus sequence drop off rapidly.

```python
def create_matrix(distances):
    matrix = np.tril(distances.values).tolist()
    for i in range(len(matrix)):
        matrix[i] = matrix[i][:i+1]
    dm = DistanceMatrix(list(distances.index), matrix)
    return dm
matrixs = [create_matrix(distances) for distances in distance_sets]
matrixs
```

### Construct a tree for each distance matrix

```python
def make_tree(dm):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree

trees = [make_tree(dm) for dm in matrixs]
trees
```

### Draw a tree for each matrix

```python
%matplotlib inline

def draw_tree(tree):
    tree.ladderize()
    Phylo.draw(tree)

for tree in trees:
    draw_tree(tree)
```

This shows us visually that branch length drops off rapidly, with only the first few longest branches being 
over 50 in length. After that point, branch length stays well below 10, and only continues to decrease.

```python

```
