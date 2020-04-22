---
permalink: /mhowland/
title: "McClane Howland"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---
[This paper](https://www.nature.com/articles/s41591-020-0820-9) explains three possible theories for the origin of the SARS-CoV-2 virus. In short, they are:

1. Natural selection in an animal host before zoonotic transfer
2. Natural selection in humans following zoonotic transfer
3. Selection during passage

Regarding hypothesis 2, the authors suggest that
"Studies of banked human samples could provide information on whether such cryptic spread has occurred."
Theories 1 and 2  also necessitate that either an animal or human virus ancestor would have a polybasic cleavage site, as does the SARS-CoV-2 virus. Existing databases of influenza viruses, [GISAID](https://www.gisaid.org/), and of coronaviruses, [CoVDB](http://covdb.popgenetics.net/v2/) could be searched for viruses that contain a polybasic cleavage site. The paper specifies the location of the cleavage site in the S protein.

# SARS-CoV-2 Spike Protein Cleavage Site

According to [this paper](https://doi.org/10.1038/s41591-020-0820-9), a polybasic cleavage site appears in the SARS-CoV-2 genome at the junction of the S1 and S2 subunits of the spike protein, which appears as amino acids RRAR at positions 682 - 685. From the [UCSC Genome browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=wuhCor1&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=NC_045512v2%3A21563%2D25384&hgsid=826092117_1DuYWre1PsFQYzKEAeA45quha3hk), the S1 protein is located at 21599 - 23617 base pairs, and the S2 protein is located at 23618 - 25381 base pairs.

This cleavage site affects the shape of the spike (is it cleaved inside the cell, outside the cell, or as it attaches to a cell?) and may be an important feature that ancestors of the virus should also have. The cleavage site also plays a role in what hosts a virus can reproduce in, because the enzymes of that host must be able to work on the cleavage site.

## Finding the cleavage site

After downloading the SARS-CoV-2 spike protein, [P0DTC2](https://covid-19.uniprot.org/uniprotkb/P0DTC2), we can see the cleavage site:


```python
from Bio import SeqIO

with open("sars-cov-2-P0DTC2.fasta") as f:
    for r in SeqIO.parse(f, "fasta"):
        print(r)
        print("Cleavage site:", r.seq[681:685])
```

    ID: sp|P0DTC2|SPIKE_SARS2
    Name: sp|P0DTC2|SPIKE_SARS2
    Description: sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=S PE=1 SV=1
    Number of features: 0
    Seq('MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDL...HYT', SingleLetterAlphabet())
    Cleavage site: RRAR


## Comparing to other proteins

I downloaded 252 BLAST alignments of the SARS-CoV-2 spike protein against various other coronaviruses and random proteins found in organisms from [uniprot.org](https://www.uniprot.org/), [blast here](https://www.uniprot.org/blast/uniprot/B20200421DA437993067D6F64326E5E763500BDED113ACDW), [alignments here](https://www.uniprot.org/align/A20200422DA437993067D6F64326E5E763500BDED01CD0BZ.aln).

Let's see if we can find one that also contains RRAR amino acids


```python
print("Using raw sequences:")
with open("B20200421DA437993067D6F64326E5E763500BDED113ACDW.fasta") as f:
    for r in SeqIO.parse(f, "fasta"):
        res = r.seq.find("RRAR")
#         print(r.description.split()[3])
        if res > -1:
            print(res, r.description)

print()
print("using alignments:")
with open("A20200422DA437993067D6F64326E5E763500BDED01CD0BZ.aln.txt") as f:
    sequences = AlignIO.read(f, "clustal")
    for s in sequences:
        r = s.seq.find("RRAR")
        if r > -1:
            print(r, s.description)
```

    Using raw sequences:
    623 sp|P11225|SPIKE_CVMJH Spike glycoprotein OS=Murine coronavirus (strain JHM) OX=11144 GN=S PE=1 SV=1
    764 sp|Q02385|SPIKE_CVMJC Spike glycoprotein OS=Murine coronavirus (strain JHMV / variant CL-2) OX=33735 GN=S PE=1 SV=1
    764 sp|P22432|SPIKE_CVM4 Spike glycoprotein OS=Murine coronavirus (strain 4) OX=12760 GN=S PE=3 SV=1
    504 tr|A0A1X9JJW0|A0A1X9JJW0_9ALPC Spike protein OS=Coronavirus AcCoV-JC34 OX=1964806 PE=4 SV=1
    40 tr|A6RAB8|A6RAB8_AJECN Uncharacterized protein OS=Ajellomyces capsulatus (strain NAm1 / WU24) OX=339724 GN=HCAG_05906 PE=4 SV=1
    
    using alignments:
    851 D6F643
    847 SP|P11225|SPIKE_CVMJH


It looks like the Murine (a rat) coronavirus spike protein also has a similar sequence of amino acids between it's S1 and S2 spike subunits. I'm not sure if the presence of an RRAR sequence directly corresponds to a cleavage site, and we can't conclude much from just this data.

## Constructing a phylogeny tree from these alignments



```python
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo


with open("A20200422DA437993067D6F64326E5E763500BDED01CD0BZ.aln.txt") as f:
    sequences = AlignIO.read(f, "clustal")
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(sequences)
    
dtconstructor = DistanceTreeConstructor()
tree = dtconstructor.nj(dm)
tree.ladderize()
Phylo.draw_ascii(tree)
```

     , SP|Q3I5J5|SPIKE_BCRP3
     |
     | SP|Q3LZX1|SPIKE_BCHK3
     |
     , TR|Q0QDX9|Q0QDX9_CVHSA
    _|
     | SP|Q0Q475|SPIKE_BC279
     |
     | __ TR|R9QTA0|R9QTA0_CVHSA
     ||
     || _ TR|A0A0K1Z074|A0A0K1Z074_CVHSA
      ||
      ||      , TR|A0A0U1UYX4|A0A0U1UYX4_CVHSA
      || _____|
       ||     | TR|A0A166ZL64|A0A166ZL64_9NIDO
       ||
       ||      , TR|D2E1D2|D2E1D2_CVHSA
        | _____|
        ||     , TR|D2E1E7|D2E1E7_CVHSA
        ||     |
        ||     | SP|P59594|SPIKE_CVHSA
         |
         | _______ D6F643
         ||
         ||           __________________ TR|A0A088DJY6|A0A088DJY6_9BETC
          |          |
          |          |            ___________ TR|A0A1B3Q5W5|A0A1B3Q5W5_9BETC
          |__________|   ________|
                     |  |        | __________ SP|A3EXG6|SPIKE_BCHK9
                     |  |        ||
                     |  |         |__________ TR|A3EXJ0|A3EXJ0_BCHK9
                     |__|
                        | ___________________ SP|P11225|SPIKE_CVMJH
                        ||
                        ||              ________ TR|A0A023Y9K3|A0A023Y9K3_9BETC
                        ||            ,|
                         |            ||       _ TR|A3EXF7|A3EXF7_BCHK5
                         |            ||______|
                         |            |       |_ SP|A3EXD0|SPIKE_BCHK5
                         |____________|
                                      |         , TR|R9UQ53|R9UQ53_9BETC
                                      |_________|
                                      |         | SP|K9N5Q8|SPIKE_CVEMC
                                      |
                                      |        , SP|Q0Q4F2|SPIKE_BC133
                                      |________|
                                               |, TR|A3EXC1|A3EXC1_BCHK4
                                               ||
                                                , SP|A3EX94|SPIKE_BCHK4
                                                |
                                                | TR|A3EXA3|A3EXA3_BCHK4
    


## Not very useful

Though the spike proteins aren't a terrible match, the whole genomes of this murine coronavirus and the SARS-CoV-2 virus are pretty different. I highly doubt the murine coronavirus is an ancestor of the current outbreak.

# Old Project 1 stuff

# [Nextstrain](https://nextstrain.org/)

Visualizations and analyses of 2019-nCoV genome evolution, using sequences from GISAID.


# [COVID-19 Galaxy Project](https://covid19.galaxyproject.org/)

The Galaxy Project is a tool for creating and publishing easily reproducible computational workflows and analyses. There is a website dedicated to sharing workflows related to the current outbreak. For example, [this](https://covid19.galaxyproject.org/cheminformatics/#virtual-screening-of-the-sars-cov-2-main-protease-de-nbi-cloud-stfc) is a workflow for screening compounds for binding to the SARS-CoV-2 protease. The SARS-CoV-2 protease is a protein that cleaves strands of amino acids formed from transcription of viral RNA into smaller pieces that will then assemble to form the virus. Finding a molecule that binds to and inhibits the action of the SARS-CoV-2 protease could help stop it's replication in cells.


# [Folding@Home](https://foldingathome.org/)

A massive distributed computing project using volunteer computers to simulate protein dynamics. Among many other ongoing projects, screenings of potential inhibitors to the same protease described above are currently being conducted. Simulations to understand how the spike protein of SARS-CoV-2 virus opens to interact with ACE2 receptors (also a protein) on human cells.

[Folding@Home's method of simulating protein folding](https://foldingathome.org/dig-deeper/#how-does-foldinghome-simulate-protein-folding)

# [DeepMind's AlphaFold](https://deepmind.com/research/publications/AlphaFold-Improved-protein-structure-prediction-using-potentials-from-deep-learning)

A project from Google that tries to predict the shape of a protein from its amino acid sequence. The method is described by Google below:

"In this work, we show that we can train a neural network to accurately predict the distances between pairs of residues in a protein which convey more about structure than contact predictions. With this information we construct a potential of mean force that can accurately describe the shape of a protein. We find that the resulting potential can be optimised by a simple gradient descent algorithm, to realise structures without the need for complex sampling procedures."
