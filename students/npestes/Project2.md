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

## Project Part 2 (a.k.a. Project 2 on Schedule)
While you are learning about evolutionary trees and distance matrices this week in lecture/lab we can get started on pulling together the data we need for our COVID-19 analysis. We can also do some starter analysis. Our goal over the next few weeks will be to clean, refine, and otherwise make detailed, explained, reproducible analysis.

The current complete genomes can be downloaded from here:
https://covid19.galaxyproject.org/genomics/4-Variation/current_complete_ncov_genomes.fasta

Disclaimers: I definitely expect for the analysis below to change. This isn't a lab. This is a first cut to get us discussing. We need to be critical of this work and look for problems and ways to improve it. We need to be skeptical scientists. We need to dig into the literature not just about the tools, but about the biology itself before putting out potentially misleading information. 

Technical Disclaimers: I've written this to run on my system. I want you to use this notebook as motivation and guidance to approach this weeks project. Some of the things I've done below may or may not be in your data analysis/programming/data science wheel house. That's ok. There is a lot of room for variety. Try to approach this in a genuine way about what you are interested in contributing.


### Guidance (i.e., what you should do)
We've said many times. This project isn't about everyone reaching the same point in a predetermined set of steps. It's about applying what we are learning in class to produce real data analysis for the community. It is about as *Learn by doing* as you could possibly get at Cal Poly. So what should you be doing this week for the project? Here is some guidance (but remember this is only to guide you and not box you into specific tasks). They are in no particular order. 
* Consider what questions we want to ask from our evolutionary tree analysis. Think about what questions the book was trying to answer. Do we even have the data in this notebook to answer some of those questions? If not, spend time trying to find it now that you can know more about what to look for in terms of format. Do some literature searching and see what other work has been done for this virus and others.
* Research and try different evolutionary tree programs/frameworks. What I've done below is not the only game in town by far. Biopython itself has different options.
* Consider the alignment itself. Are there different ways to do this? Did we do it correctly?
* What about the sequences themselves? Are they all of the same quality? Should we exclude some?
* What about the virus alignment program? Did we use that correctly? Should we have done the entire sequence instead of using Spike as a reference? Should we try a different reference. 
* Do we have more data available about the sequences? Part of world, etc. Can we do some digging here to answer different questions.
* And I'm sure you can think of more to attempt... Think about what you want to do. Spend time working towards a well thoughtout goal. Document things as you go. Talk to everyone on Slack. Together we can do this!


### Link to clone the repository
Here is a link to the project repository.

https://github.com/anderson-github-classroom/csc-448-project

The website can be viewed at https://anderson-github-classroom.github.io/csc-448-project/.


### First step is to get the data
We are going to rely on the Galaxy team to pull together our sequence data for now. We might change this later.

```python
import wget

url = 'https://covid19.galaxyproject.org/genomics/4-Variation/current_complete_ncov_genomes.fasta'
file = '../../current_complete_ncov_genomes.fasta'
wget.download(url, file)
```

### Virus Alignment
We will use a program specific for viral multiple alignments: https://github.com/rega-cev/virulign-tutorial

https://academic.oup.com/bioinformatics/article/35/10/1763/5123354

I downloaded the Mac binary and put it /Users/panderson/

```python
!git clone https://github.com/rega-cev/virulign.git
```

```python
!ls
!ls virulign/build/src/Debug/virulign.exe
```

I also downloaded the tutorials and the program repository.

```python
!git clone https://github.com/rega-cev/virulign-tutorial ../../virulign-tutorial
```

```python
!git clone https://github.com/rega-cev/virulign ../../virulign
```

### Before alignment
As we mentioned in class, we need an alignment so we can derive our pairwise distance scores so we can then put together our distance matrix.

This package contains a reference Spike protein that can be provided as an argument when performing alignment. This code took my computer a few minutes to run, so I've included the output in the project repository: csc-448-project/data/position_table.csv.

```python
!virulign/build/src/Debug/virulign.exe ../../virulign/references/SARS-CoV-2/S.xml ../../current_complete_ncov_genomes.fasta --exportAlphabet Nucleotides --exportKind PositionTable > ../../position_table.csv
```

### Read the data into a pandas dataframe

```python
import pandas as pd
position_table = pd.read_csv('../../data/position_table.csv') # or put in the path to csc-448-project/data/position_table.csv
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
import numpy as np
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


Try with all sequences

```python
subset_seqs2 = distance_from_concensus_seq_sorted.index
subset_seqs2
```

```python
distances = {}
for i,seqid1 in enumerate(subset_seqs2):
    distances[seqid1,seqid1]=0
    for j in range(i+1,len(subset_seqs2)):
        seqid2 = subset_seqs2[j]
        distances[seqid1,seqid2] = sum(position_table.loc[seqid1] != position_table.loc[seqid2])
        distances[seqid2,seqid1] = distances[seqid1,seqid2]
distances = pd.Series(distances).unstack()
distances
```

```python
from Bio.Phylo.TreeConstruction import DistanceMatrix
import numpy as np
matrix = np.tril(distances.values).tolist()
for i in range(len(matrix)):
    matrix[i] = matrix[i][:i+1]
dm = DistanceMatrix(list(distances.index), matrix)
```

```python
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
```

```python
%matplotlib inline

from Bio import Phylo
tree.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(tree)
```

Lots of very similar strains near the base


try with UPGMA instead of NJ

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

```python
from Bio.Phylo.TreeConstruction import DistanceMatrix
import numpy as np
matrix = np.tril(distances.values).tolist()
for i in range(len(matrix)):
    matrix[i] = matrix[i][:i+1]
dm = DistanceMatrix(list(distances.index), matrix)
```

```python
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
```

```python
from Bio import Phylo
tree.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(tree)
```

Distinctly different tree generated fewer strains shown close to the base of the tree

```python
subset_seqs3 = distance_from_concensus_seq_sorted[:40].index
subset_seqs3
```

```python
distances = {}
for i,seqid1 in enumerate(subset_seqs3):
    distances[seqid1,seqid1]=0
    for j in range(i+1,len(subset_seqs3)):
        seqid2 = subset_seqs3[j]
        distances[seqid1,seqid2] = sum(position_table.loc[seqid1] != position_table.loc[seqid2])
        distances[seqid2,seqid1] = distances[seqid1,seqid2]
distances = pd.Series(distances).unstack()
distances
```

```python
from Bio.Phylo.TreeConstruction import DistanceMatrix
import numpy as np
matrix = np.tril(distances.values).tolist()
for i in range(len(matrix)):
    matrix[i] = matrix[i][:i+1]
dm = DistanceMatrix(list(distances.index), matrix)
```

```python
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()
treenj = constructor.nj(dm)
treeupgma = constructor.upgma(dm)
tree.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(treenj)
```

```python
tree.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(treeupgma)
```

MT233522.1 appears as the farthest out in nj, but as the earliest in upgma


We have examined variations in the spike protein of COVID19. We have used distance metrics that treat a changes in any position of the protein as equally important. While reasonable for phylogenetic research it is not true that all changes are equally likely, or that all changes are equally important.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7151553/ provides a recent summary of mutation patterns in the COVID19 spike protein. Targeting conserved regions may prove an effective therapeutic, since the slowed mutation may indicate important (difficult to change) sections of the protein.

```python
import pandas as pd
position_table = pd.read_csv('../../data/position_table.csv') # or put in the path to csc-448-project/data/position_table.csv
```

```python
num_mutations = position_table.nunique()[1:]
```

```python
max(list(num_mutations))
```

```python
%matplotlib inline
import matplotlib.pyplot as plt

plt.hist(num_mutations)
```

```python
total = 7+356+3456
print(7/356)
print(356/3456)
print(7/total, 356/total, 3456/total)
```

```python
expected3 = (356/total)*(356/total)*total
expected2 = (356/total)*total
expected1 = (1-356/total)*total
```

```python
xi2 = ((7 - expected3)**2)/expected3 + ((356-expected2)**2)/expected2 + ((3456-expected1)**2)/expected1
xi2
```

```python
(((7 - expected3)**2)/expected3, ((356-expected2)**2)/expected2, ((3456-expected1)**2)/expected1)
```

the xi squared statistic for 1% chance of occurance with degree of freedom 2 (3-1) is 9.21. We exceed this. So we reject the null with greater than 99% confidence.
We used the null hypothesis of a exponential distribution that perfectly matches the distribution at 2

```python
import math

prob = math.sqrt((7/total))
expected3 = prob*prob*total
expected2 = prob*total
expected1 = (1-prob)*total
```

```python
xi2 = ((7 - expected3)**2)/expected3 + ((356-expected2)**2)/expected2 + ((3456-expected1)**2)/expected1
xi2
```

```python
prob = 1-(3456/total)
expected3 = prob*prob*total
expected2 = prob*total
expected1 = (1-prob)*total
```

```python
xi2 = ((7 - expected3)**2)/expected3 + ((356-expected2)**2)/expected2 + ((3456-expected1)**2)/expected1
xi2
```

```python
other reasonable hypotheses are equally bad
```

```python
import numpy as np
gaus_size = 9
alpha = 0.1
gaus_filter = np.exp(-alpha*(np.arange(gaus_size) - gaus_size/2)**2)
gaus_filter = gaus_filter / np.sum(gaus_filter)
gaus_filter = gaus_filter / np.sum(gaus_filter)
local_mutations = np.convolve(np.array(list(num_mutations)), gaus_filter)
plt.hist(local_mutations)
```

```python
gaus_filter
```

```python

```
