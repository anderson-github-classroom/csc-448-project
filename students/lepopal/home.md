---
permalink: /lepopal/
title: "Lemar Popal"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
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


Here is a notebook where I download all the Covid-19 sequence files from the NCBI database: [Getting Data](GettingData.md)

## Week 2
