---
permalink: /lepopal/
title: "Lemar Popal"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---

# Lemar Popal

## Week 1
First conducting some research...

https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/

Here’s a database of nucleotide sequences of COVID-19 collected at different dates throughout the world. I think this might be useful to see how the virus has mutated over time, and possibly see which or how many variations of the virus there were at any given time. 

---

[biopython.org](biopython.org)

http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc117

This looks like one of the best Python libraries out there for biology. It seems like it includes a lot of useful tools, but most of them I don’t know to use (yet!). Using it I was able to relatively easily download all genome sequences from the NCBI database into Sequence objects in Biopython. From there we can do a lot of things, like finding the reverse complement using the built in method in the Sequence object. 

---

TODO: Conduct phylogenetic tree, don't have to analyze yet 

```python
import pandas as pd
import Bio as bp
import numpy as np
import re
```

First, I want to download the various (complete) Coronavirus genomes. It might be useful later.

From the NCBI database I downloaded a .yaml file which contained the "accession" ID of every (complete) sequence. Using these IDs I can download each sequence using Biopython. I have to do some cleanup of the file so I can get the IDs I want.

```python
ids = []
with open("ncov-sequences.yaml","r") as f:
    lines = [line.strip() for line in f.readlines() if line not in 'genbank-sequences:\n'][1:]
```

```python
# get indices of all lines with 'accession' so we know where to slice list
indices = [i for i, line in enumerate(lines) if '- accession' in line]
```

```python
# get all sequences
sequences = []
for i in range(len(indices)):
    try:
        sequences.append(lines[indices[i]:indices[i+1]])
    except IndexError:
        sequences.append(lines[indices[i]:indices[i]+6])
```

```python
# filter for sequences that are complete
complete_sequences = []
for lst in sequences:
    for line in lst:
        if "gene-region: complete" in line:
            complete_sequences.append(lst)
```

```python
# use regex to grab the ids of the complete sequences
ids = []
for i in range(len(indices)):
    match = re.search('accession: [A-Z]{2}[\d]{6}',lines[indices[i]])
    if match:
        ids.append(match.group(0)[-8:])
```

The above code was probably pretty ugly, but now we have a list of ids that we can query the database with. 

```python
import os
from Bio import SeqIO
from Bio import Entrez

for _id in ids: 
    Entrez.email = "lepopal@calpoly.edu"  
    filename = str(_id)+".gbk"
    if not os.path.isfile("sequences/" + filename):
        # Downloading...
        net_handle = Entrez.efetch(
            db="nucleotide", id=_id, rettype="gb", retmode="text"
        )
        out_handle = open("sequences/" + filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved", filename)
    else:
        print("Already exists", filename)
```

Now we have all the files. We can load them all into Biopython. 

```python

```
