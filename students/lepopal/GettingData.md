---
jupyter:
  jupytext:
    formats: ipynb,md
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

# Lemar Popal

## Week 1
My goal this week was to do some research on the various resources already out there on/about COVID-19 (the coronavirus). By researching what people have done already I can hopefully come up with ideas I can research and implement myself and be able to contribute meaningfully to the scientific community. 

---

[https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/)

Here’s a database of nucleotide sequences of COVID-19 collected at different dates throughout the world. I think this might be useful to see how the virus has mutated over time, and possibly see which or how many variations of the virus there were at any given time. 

---

[biopython.org](biopython.org)

This looks like one of the best Python libraries out there for biology. It seems like it includes a lot of useful tools, but most of them I don’t know to use (yet!). Using it I was able to relatively easily download all genome sequences from the NCBI database into Sequence objects in Biopython. From there we can do a lot of things, like finding the reverse complement using the built in method in the Sequence object. 

---

[gisaid.org](https://www.gisaid.org/)

This site also looks like it will be useful in my research. They make it easy for scientists to share information about virus sequences, clinical/epidemiological data, and geographic data associated with viruses. I found that [nextstrain.org/ncov/](https://nextstrain.org/ncov/) uses data from Gisaid. The visualizations include the phylogenic tree of the coronavirus (rooted to early samples of the virus from Wuhan) and shows a map of how the virus was spread and how the virus mutated over ttime. I think once I learn more about bioinformatics in the coming weeks I'll be able to conduct some new analysis that this website hasn't done yet. 

---
[https://covid19.galaxyproject.org/](https://covid19.galaxyproject.org/)

This goal of this site is to "provide publicly accessible infrastructure and workflows for SARS-CoV-2 data analyses." Specifically, there are three different types of analysies that they feature: genomics, evolution, and cheminformatics. A lot of these methods I don't know how to use yet, but I anticipate that we'll touch on them later in the book and that I'll be able to access these workflows, implement them, and come up with some meaningful results. 

---

TODO: Create a phylogenetic tree, don't have to analyze yet.


First, I want to download the various (complete) Coronavirus genomes. It might be useful later.

From the NCBI database I downloaded a .yaml file (available [here](https://www.ncbi.nlm.nih.gov/core/assets/genbank/files/ncov-sequences.yaml)) which contained the "accession" ID of every sequence in their database. Using these IDs I can download each sequence using Biopython. I have to do some cleanup of the file so I can get the IDs I want. 

```python
import Bio as bp
```

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
import re

# use regex to grab the ids of the complete sequences
ids = []
for i in range(len(indices)):
    match = re.search('accession: [A-Z]{2}[\d]{6}',lines[indices[i]])
    if match:
        ids.append(match.group(0)[-8:])
```

The above code was probably pretty ugly, but now we have a list of IDs that we can query the database with. 

```python
import os
from Bio import SeqIO
from Bio import Entrez

for _id in ids[:10]: # remove the slice notation to get all IDs. This will only grab first 10 sequences. 
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
        print("Saved:", filename)
    else:
        print("Already exists:", filename)
```

Now we have all the files. We can load them all into Biopython. 

```python
files = [file for file in os.listdir("sequences/") if file.endswith(".gbk")]

for fname in files:
    for seq_record in SeqIO.parse("sequences/" + fname, "genbank"):
        print(seq_record.id) # prints the ID of the sequence
        print(repr(seq_record.seq)) # prints the nucleotide sequence
        print(len(seq_record)) # prints the length of the nucleotide sequence
```

At this point I think I have a good structure to start working with the data in these files. Using what we learned in Week 1, next week I want to try to find the *ori* of the virus. Also, I want to try and create a phylogenetic tree like on [nextstrain.org/ncov/](https://nextstrain.org/ncov/). 


## Week 2

```python

```
