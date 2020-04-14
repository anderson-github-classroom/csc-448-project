---
permalink: /lepopal/
title: "Lemar Popal"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---

### Week 1
https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/

Here’s a database of nucleotide sequences of COVID-19 collected at different dates throughout the world. I think this might be useful to see how the virus has mutated over time, and possibly see which or how many variations of the virus there were at any given time. 

---

[biopython.org](biopython.org)

This looks like one of the best Python libraries out there for biology. It seems like it includes a lot of useful tools, but most of them I don’t know to use (yet!). Using it I was able to relatively easily download all genome sequences from the NCBI database into Sequence objects in Biopython. From there we can do a lot of things, like finding the reverse complement using the built in method in the Sequence object. 

```python
import pandas as pd
import Bio as bp
```

First, I want to download the various (complete) Coronavirus genomes. It might be useful later. 

```python
seq_record = bp.SeqIO.parse("sequence.fasta", "fasta")
print(seq_record.seq)
```



```python

```
