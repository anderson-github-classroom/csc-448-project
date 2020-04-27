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

## Project Part 2 

This Part of the Project will use the complete COVID-19 Genome. This genome was provided by the galaxy project, and can be found here: https://covid19.galaxyproject.org/genomics/4-Variation/current_complete_ncov_genomes.fasta


### Tools
Biopython
Biopython was utilized during this assignment, and was installed using
pip install biopython

Also used in this project was virulign. Virulign can be found at: https://github.com/rega-cev/virulign


### Pulling Galaxy Project Data
The previously cited Galaxy project data was gathered using the following code:

```python
import wget

url = 'https://covid19.galaxyproject.org/genomics/4-Variation/current_complete_ncov_genomes.fasta'
file = '../../current_complete_ncov_genomes.fasta'
wget.download(url, file)
```

### Virus Alignment
When running virulign, I decided to compare the genome to the HIV-HXB2-env genome provided by virulign, because of all the rumors that have been spreading about the disease having a similar structure to HIV. I wanted to test how close the relation may be, and used the following commands to do so.

.\virulign .\HIV\HIV-HXB2-env.xml current_complete_ncov_genomes.fasta > alignment.mutations 2> alignment.err


### Read the data into a pandas dataframe

```python
import pandas as pd
position_table = pd.read_csv('../../position_table.csv') # or put in the path to csc-448-project/data/position_table.csv
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

### Select 10 sequences to do our first analysis

```python
subset_seqs = distance_from_concensus_seq_sorted[:10].index
subset_seqs
```

### Construct a distance matrix for our sequences

```python
distances = {}
for i,seqid1 in enumerate(subset_seqs):
    distances[seqid1,seqid1]=0
    for j in range(i+1,len(subset_seqs)):
        seqid2 = subset_seqs[j]
        distances[seqid1,seqid2] = sum(position_table.loc[seqid1] != position_table.loc[seqid2])
        distances[seqid2,seqid1] = distances[seqid1,seqid2]
distances = pd.Series(distances).unstack()
distances
```

### Utilize biopython
For this analysis we'll use a package called biopython: ``pip install biopython``. 

It has its own formats, so we'll need to convert.

```python
from Bio.Phylo.TreeConstruction import DistanceMatrix
matrix = np.tril(distances.values).tolist()
for i in range(len(matrix)):
    matrix[i] = matrix[i][:i+1]
dm = DistanceMatrix(list(distances.index), matrix)
```

### Now construct our tree

```python
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
```

### Now draw our tree

```python
%matplotlib inline

from Bio import Phylo
tree.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(tree)
```

**Please see the guidance at the top of the page for what to try**

```python

```
