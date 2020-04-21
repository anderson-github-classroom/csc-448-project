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
> Studies of banked human samples could provide information on whether such cryptic spread has occurred.
Theories 1 and 2  also necessitate that either an animal or human virus ancestor would have a polybasic cleavage site, as does the SARS-CoV-2 virus. Existing databases of influenza viruses, [GISAID](https://www.gisaid.org/), and of coronaviruses, [CoVDB](http://covdb.popgenetics.net/v2/) could be searched for viruses that contain a polybasic cleavage site. The paper specifies the location of the cleavage site in the S protein.

Platform dedicated to sharing flu virus genomes. For example, [this paper](https://mra.asm.org/content/9/11/e00169-20) used a reference sequence from GISAID and uploaded their sequence to GISAID after.

# Old stuff

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

