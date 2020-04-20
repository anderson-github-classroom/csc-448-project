---
permalink: /wiobrien/
title: "William O'Brien"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---

# Week 1
## The GISAID Initiative
### https://www.gisaid.org/about-us/mission/

GISAID is a non-profit initiative recognized by the European Commission as a parter in the PREDEMICS (Prepardness, Prediction and the Prevention of Emerging Zoonotic Viruses with Pandemic Potential) consortium. As such, their research regarding COVID19 has been multinational and exhaustive. A number of the following sources utilize GISAID data in their own analysis, emphasizing its importance in the field of bio-informatics. At the time of writing, they have a total of 8,001 genomes sequenced regarding hCoV-19. Unlike conventional archives, after curation, genome sequences are released to the public immediately. 

## Folding @ Home Project
### https://foldingathome.org/

The Folding @ Home Project is based in St. Louis at the Washington University in St. Louis School of Medicine. Their research investigates the protein folds in human cells which are integral to ensuring proper cell functionality. Misfolds in these proteins can result in serious alterations of functionality, possibly introducing serious health complications. Regaring COVID19, Folding @ Home is using computer simulations to understand the viral proteins found in instances of COVID19. Users can donate their own computing power to the research of these folds by visiting the URL above. 

## Next Strain
### https://nextstrain.org/ncov

Next Strain aggregates data into a tree-based representation of the fluctuation of COVID19 as it spreads into new territories. The distance of a particular node from the root correlates with the number of mutations of that particular instance of COVID19. The default view shows 3414 genomes but far more can be viewed if the user narrows the time/geographical ranges. The parent site, nextstrain.org, contains a plethora of data and analysis from a broad range of international sources. Analyses and reports can be found for each datasource and genome.

## Complete Genome Sequence of a 2019 Novel Coronavirus (SARS-CoV-2) Strain Isolated in Nepal
### https://mra.asm.org/content/9/11/e00169-20

A complete genome sequence of an instance of COVID-19 can be viewed from the data embedded in this study. This particular strand originates from a 32-year old patient who contracted COVID-19 in Wuhan, China. He then travelled to Nepal where he was diagnosed and studied. Research on this full genome can help research the evolution of COVID-19 and aid in the development of vaccines and treatment. The data is readibly available for individual research.

## Best practices for the analysis of SARS-CoV-2 data: Genomics, Evolution, and Cheminformatics
### https://covid19.galaxyproject.org/

The Galaxy Project analyzes data from GISAID in order to investigate the genomics, evolution, and cheminformatics of COVID-19. The site uses open-source data-science tools to provide easily-modifiable templates and workflows for analysis on each of the fields above. These templates can be downloaded and modified easily, making best-practice tools available for public use.

# Week 2
## Relative Distance Analysis to Suggest Host of Transfer for COVID19

The successing cells load data from the NCBI Nucleotide database as a means of investigating the likely host that humans initally obtained COVID19 from. With a relatively small sample size of host species we can identify how significant the edit distance between the sequenced human genome and the sequenced bat genome, further bolstering the current theory that COVID19 was originally transferred from bats to humans.

### Data Collection
I used a set amount of accessions from a varying source of hosts for my data set. The hosts were arbitrarily chosen from the NCBI nucleotide database.

```python
accessions = ["1798174254", "MN996532","NC_002306","KY566211" 
              ,"JF792617", "KM454473", "NC_017083", "KR061459","KF294357", 
              "NC_003045", "LN610099", "MH021175", "KF268339",
             "NC_003436"]
host_names = ["human", "bat", "feline", "feline2", 
              "rat", "duck", "rabbit", "swine", "mouse", 
              "bovine", "guinea fowl", "avial (Brazil)", "murine",
             "porcine"]
```

### Data Source
I used the Entrez subset of the Bio module for data collection. I used the nucleotide DB within the NCBI database list to obtain full nucleotide sequences of the above accessions.

```python
from Bio import Entrez
from Bio import SeqIO
handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="xml")

record = Entrez.read(handle)

handle.close()

```


This data was then placed into a blank pandas DataFrame in a 2D matrix.

```python
import pandas as pd
place_holder = []
for i in range (len(host_names)):
    place_holder.append([0 for j in range(len(host_names))])
print(place_holder)
D = pd.DataFrame(place_holder,index=host_names,columns=host_names)

D
```
![image](dm.png)

The distances between nucleotides were computed using the python Levenshtein module. The standard Levenshtein distance algorithm was employed in this case. For large accession counts, this process can take a few minutes.

```python
from Levenshtein import *

for name_a in host_names:
    for name_b in host_names:  
        if name_a == name_b:
            continue
        if D[name_a][name_b] != 0 or D[name_b][name_a]:
            continue
        dist = distance(sequences[name_a], sequences[name_b])
        D.set_value(name_b,name_a, dist)



# Create a lower-triangular matrix
j = 0
as_num = D.to_numpy().tolist()
tri_dm = []
for j in range(len(as_num)):
    tri_dm.append(as_num[j][0: j + 1])
print(tri_dm)

D
```

Finally, I displayed the graph to showcase the evolutionary tree of the selected coronavirus sequences.

```python
%matplotlib inline

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
constructor = DistanceTreeConstructor()
tree = constructor.nj(DistanceMatrix(names=host_names, matrix=tri_dm))
tree.ladderize() 
Phylo.draw(tree)
```

![image](wk2tree.png)

### Analysis
As you can see, the human and bat genomes feature a statistically significant lack of variation relative to the other coronavirus sequences selected. This claim is bolstered by the distance matrix as well which showcases just how similar the two host sequences are. The current theory on the transmissin of COVID19 suggests that the virus was first contracted by humans through contact with bats. This very simple demonstration reaches a similar conclusion, albeit with a relatively small sample size.

### Room for Improvement
There is ample room for improvement in this particular analysis. For starters, relatively few accessions were selected. A better approach would be to take a large, random sample from the set of hosts in the NCBI database that have entire coronavirus genomes sequenced. Additionally, sequencing and analyzing the entire genome may take into account insignificant variation that could otherwise be filtered by analyzing particular proteins. This introduces additional challenges, however, such as ensuring that each host has the same protein sequenced. Either way, this basic research covers a lot of the tools for analysis and data collection that will prove useful in the further analysis of COVID19 and future coronaviruses.