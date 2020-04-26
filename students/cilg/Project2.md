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
!/Users/panderson/bin/virulign
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

```python jupyter={"outputs_hidden": true}
!/Users/panderson/bin/virulign ../../virulign/references/SARS-CoV-2/S.xml ../../current_complete_ncov_genomes.fasta --exportAlphabet Nucleotides --exportKind PositionTable > ../../position_table.csv
```

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
