---
permalink: /bfkwong/
title: "Bryan Kwong"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---

## I. PROJECT 1 

<!-- #region -->
This week, I will explore the available sources of bioinformatics data that are currently available on the internet right now. Since COVID-19 was first discovered in Hubei at the end of 2019, labs all across the world has been working hard to sequence the genome of this Coronavirus. While this is hard and expensive work, it is crucial for tracking how the disease is spreading across the globe as well as identify methods of intervention.

Since much work has already been done in identifying and analyzing COVID-19's genome. It is worth the effort to read through and catch up on what has already been done and build on top of that work. Fortunately in the age of the internet, a lot of this research is already online, making it easy find other research. 

### A. GISAID Initiative

The GISAID Initiative stands for Global Initiative on Sharing All Influenza Data. It was set up by the German Government and nonprofits to provide public access to the most complete genome data on influenza and other diseases (such as COVID-19). GISAID provides valuable tools such as tracking the epidemiology of Influenza to try to predict which strain is most likely to occur as well as provide recommendation for which strain of Influenza we should vaccinate against in the up coming flu season. 

These tools were incredible valuable for influenza but due to the modular nature of GISAID, it is currently being applied to COVID-19. GISAID provides a fantastic chart about the epidemiology of COVID-19 showing how mutations occured from the original strain discovered in late 2019. Also, they provide a dashboard that shows confirmed cases, active cases, and the outcome of resolved cases. From GISAID and the information they provide, I learned the following: 

**CASES OF CORONAVIRUS**: As of 4/7/2020, there are 1.4 million confirmed cases of COVID-19 with 81,000 deaths and 300,000 recovered. Most of the new cases of COVID-19 are outside of China (Europe and the US). Currently, the US has the most confirmed cases of COVID-19 at 398,000 cases and Italy with the most deaths at 17,000 deaths.

**EPIDEMIOLOGY OF COVID-19 IN THE US**: GIDAID's epidemiology chart allows us to track and see when and where the viruses were identified in the US. Additionally, it also tells us where it most likely came from. The following are some key events when it comes to COVID-19 arriving in the US. 

1. **MN985325**: Man travelled from China diagnosed in Snohomish County, WA
2. **MN988713**: Man travelled from China diagnosed in Chicago, IL 
3. **MN994468**: Man travelled from China diagnosed in Orange County, CA
4. **Mutation occurred in Europe**: Major mutation occured in Europe on 1/16/2020
5. **419553**: First person diagnosed with COVID-19 that traced back to Europe in Rhode Island
6. **Currently two clusters of COVID-19 Phylogeny in the US.** One cluster originating from China, the other originating from Europe

### B. Folding at Home

Folding at Home is a project sponsored by Oracle that allows people to pool together researches for protein folding simulations. These simulations are useful because they allow us to see how viruses behave and to observe protein behaviors that often goes undetected in lab experiments. An example they gave was Ebola, where the team was able to find an opening that could potentially allow therapeutics to take effect. Folding at Home hopes to accomplish the same for COVID-19 where they can discover therapeutics for the Coronavirus through protein interactions otherwise unobservable in laboratory experiments.

### C. Complete Genome Sequence of a SARS-CoV-2 Strain Isolated in Nepal

This paper (published by the American Society of Microbiology) explained the process undertaken by the University of Hong Kong to sequence the Genome of COVID-19 (then still known as SARS-CoV-2). The Chemical process of sequencing the genome (rRT-PCR developed at the University of Hong Kongs) is complicated process that is quite difficult to understand, however the results are quite simple to interpret.  

1. **SARS-CoV-2 Nucleotide Sequence Length:** 29,811 nucleotides
2. **Nucleotide Composition:** 
  * 29.86% adenosines
  * 18.39% cytosines
  * 19.63% guanines
  * 32.12% thymines. 
  
  
3. **Mutations from COVID-19 Sequenced in Wuhan**: 
  * ORF1a, codons AGT to AGC, silent mutation
  * ORF1a, codons TTA to TCA, nonsilent mutation
  * ORF1b, codons CTA to TTA, silent mutation
  * ORF8b, codons TCA to TTA, nonsilent mutation
  * nucleocapsid, codons TTT to TTC, silent mutation
  
  
4. **Deamination Location from COVID-19 Sequenced in Wuhan**: 
  * Site 24019 vs 2019-nCoV WHU01
  * Site 24019 vs Wuhan-Hu-1
  
  
5. **Similarity with COVID-19 Sequenced in Wuhan**: > 99.99%

### D. Conclusion 

Since December 2019, a huge amount of research has been done on this Coronavirus making it one of the most thoroughly studied virus. However, there is still a lot of work to be done to fully understand how the virus works and what we can do to stop it. A lot of information is floating out there on the internet yet to be analyzed and understood. It is the purpose of CSC448 to both teach me how to do that type of analysis as well as contribute to the scientific community.

---

## I. PROJECT 2 

When a new disease emerges, one of the first things we do is to identify where the disease came from and when the disease crossed over to humans. With the 2002-2003 SARS outbreak, scientists used a phylogenetic tree to trace where the virus came from. Using this technique, bioinformaticians hypothesized that the SARS virus likely came from civet cats who likely got it from a bat. Below, I will run a similar analysis but on SARS-CoV-2.

### A. Data Retrieval

I got all my data from the NCBI Nucleotide Database queried with the Entrez tool by their accession number. I looked up their accession number by searching for Coronaviruses at the Nucleotide database on their web interface https://www.ncbi.nlm.nih.gov/nuccore. I searched by different types of viruses and the following are the viruses I choose with their accession number and host animal. 

index	|accessions|	host	|sequence
---|---|---|---|
0	|NC_009019.1|	Tylonycteris	|GATTTAAGTGAATAGCCTAGCTATCTCACCCCCTCTCGTTCTCTTG...
1	|GQ477367.1	|canine	|ACTTTTAAAGATAAAGTGAGTGTAGCGTGGCTATCTCTCATCTTTT...
2	|MK994937.1	|swine	|ACATGGGGACTTAAAGATATAATCTATCTGCCGATAGAGTCCTTAT...
3	|KX722530.1	|feline	|GTGAGTGTAGCGTGGCTATAACTCTTCTTTTACTTTAACTAGCTTT...
4	|U00735.2	|bovine	|GATTGTGAGCGATTTGCGTGCGTGCATCCCGCTTCACTGATCTCTT...
5	|KF268339.1|	murine	|TACCCTCTCAACTCTAAAACTCTTGTAGTTTAAATCTAATCTAAAC...
6|NC_009657.1	|Scotophilus|	GACTTAAAGATATTATCTATCTATAGATAGATCAATTTCTTTCCTA...
7|	JQ989273.1	|Hipposideros|	ACATGGGGACTTATTGTGATTTTCTATCTGCGGATAGTAAGTGTCA...
8|	NC_017083.1	|rabbit	|GATTCTGAGCGATTTGCGTGCGTGCATCCGCCTCAGTGAACTCTTG...
9	|MK359255.1	|goose	|ACTTTGAGCATTGATATATATATATATATATCATACTCACCTTGCC...
10|NC_032730.1	|rat	|ACTTTTAGAGTATAATCTATTATACATAGATTTGCACTAACCCCTC...
11|	NC_010438.1	|Miniopterus|	GACTTAAAGATATAATCCATCTAGCGACGGGTTATACTCTTTTTAG...
12|	LC215871.1	|ferret	|AGTGAGTGTAGCATAGCTGCCTACTTTCTTTAACTTGACTCTAAGT...
13	|NC_023760.1|	mink|	ACTTTTAAAGATAAGTGCGTGTAGCGTAGCTGCCCACCTTTCTTTA...
14	|EF424621.1|	antelope|	GCGTGCGTGCATCCCGCTTCACTGATCTCTTGTTAGATCTTTTCAT...
15	|GQ427176.1|	turkey	|ACTTAAGATAGATATTAATATATATCTATCAAACTAGCCTTGCGCT...
16	|MT072864.1	|pangolin	|TCCCAGGTAGCAAAACCAACCAACTCTCGATCTCTTGTAGATCTGT...
17	|AY572034.1	|civet	|AAGCCAACCAACCTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC...
18	|KF793825.1	|dolphin|	CTATTTGTGAATAATATATATATATATATCATTCATTCGTTTACCC...
19|	MN996532.1	|bat	|CTTTCCAGGTAACAAACCAACGAACTCTCGATCTCTTGTAGATCTG...
20|	MG518518.1	|waterdeer|	GATTGTGAGCGATTTGCGTGCGTGCATCCCGCTTCACTGATCTCTT...
21	|KM454473.1	|duck	|TGCTGGTATCACTGCTTGTTAGGTTGTGTCTCACTTTATACATCCG...
22	|NC_045512.2|	covid19	|ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGA...
23	|NC_009020.1|	Pipistrellus|	GATTTAAGAGAATAGCCTAGCTATCCCTCTCTCTCGTTCTCTTGCA...
24|	MN996532.1	|Rhinolophidae	|CTTTCCAGGTAACAAACCAACGAACTCTCGATCTCTTGTAGATCTG...
25|	NC_018871.1	|Rousettus	|ACATGGGGACTTATTGTGATTTTCTATCTGCGGATAGTAAGTGTCA...

I will use these nucleotide sequences to create a phylogenetic tree to try to determine/re-confirm which species the virus originated from. 

Instead of just using the Spike Protein Sequence like the lab, I will use the entire Nucleotide sequence which averages a length of 29472 nucleotides long. Although the textbook warned that using the entire sequence could be tricky as there are often rearrangements, insertions and deletions. However, complete sequences are the easiest set of data to compile. 

### B. Distance Matrix

The distance matrix was constructed using Levenshtein Distance. It is rather a length exercise derived from CSC349 so instead of writing it, I used the python-levenshtein package. I considered using the Damerau-Levenshtein distance but transposition added to a longer compute time, thus I settled for the Levenshtein distance. Levenshtein is the most optimal choice since or similarity measure in lab doesn't apply for my dataset. This dataset are all sequences of different length and thus Levenshtein distance was necessary. 

Since the sequence is the entire genome, there are a lot of noise in the data. Some of that noise could come from random mutations and just differences in the coronavirus type. It is my hope that by using Levenshtein distance, it could account and adjust for the random mutation so that it still gives us good results at the end. The code I used to generate the distance matrix is as follows: 

```
d_mtrx = np.zeros((len(align_str),len(align_str)))
for i in range(len(align_str)): 
    for j in range(i+1, len(align_str),1):
        dld = Levenshtein.distance(align_str.iloc[i], align_str.iloc[j])
        d_mtrx[i][j] = dld
        d_mtrx[j][i] = dld
```

### C. Phylogenetic Tree

I used BioPython to create the Phylogenetic tree using the DistanceTreeConstructor. The code for it is as follows: 

```
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
tree.ladderize() 
Phylo.draw(tree)
```

The following three phylogenetic tree was created: 

**Tree with Non-Bat Species**

![Non-Bat species](https://raw.githubusercontent.com/anderson-github-classroom/csc-448-project/master/students/bfkwong/half.png)

**Tree with Just Bat Species**

![bat species](https://raw.githubusercontent.com/anderson-github-classroom/csc-448-project/master/students/bfkwong/bat.png)

**Tree with All Species**

![all species](https://raw.githubusercontent.com/anderson-github-classroom/csc-448-project/master/students/bfkwong/full.png) 

### D. Conclusions 

From all these Phylogenetic chart that was generated. We learned / confirmed that COVID-19 likely originated from bats, specifically the horseshoe bats. A study that came out in January 27th of 2020 also stated that the likely origins of COVID-19 is from bats. 

This phylogenetic chart also showed me something very similar, which is that SARS and COVID-19 shared very similar origins. Both from horseshoe bats along the same evolutionary branch. This explains so much of why SARS and COVID-19 share such similar symptoms, and how it disproportionately affects different age groups. 
<!-- #endregion -->

```python

```
