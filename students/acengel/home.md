---
permalink: /acengel/
title: "Alexander Engel"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---

## Week 1

In order to prepare for some of the COVID-19 data we'll be looking at, let's investigate some sites.

### [GISAID Initiative](https://www.gisaid.org/)

It looks like the GISAID Initiative jumps through some hurdles to make it easy to share data on viruses.  They focus on influenza, but have also put together a collection of COVID-19 data.  Their site contains tons of different genome sequences, as well as a phylogeny and geographic map showing how the virus has mutated and spread.

This site will be very useful as a source of data.

### [Folding@Home](https://foldingathome.org/)

Folding@Home is a neat distributed computing project.  Once proteins are synthesized in our cells, they are then folded into the proper structures to be useful.  Simulating protein folding is computer-intensive, so the Folding@Home project allows volunteers to run these simulations on their idle personal devices and contribute to research.

Viruses also produce proteins that weaken our immune systems and replicate the virus.  Folding@Home is working to understand how the coronavirus proteins work in order to design therapeutics to stop them.

I'm not sure if I'll participate in Folding@Home, but I'll be sure to check out their research insights.

### [Next Strain](https://nextstrain.org/ncov)

Nextstrain is an open-source project that takes pathogen genome data and creates visualization tools for the community.  They put together a phylogeny showing how coronavirus has mutated, as well as a map showing how it has spread.  There are a ton of options for displaying the data.  They get the data from GISAID, and GISAID links to Nextstrain as an example of what can be done with their data.

This site will be useful for exploring existing visualization tools for COVID-19.

### [ASM](https://mra.asm.org/content/9/11/e00169-20)

The American Society for Microbiology (ASM) published a complete genome sequence of COVID-19 taken from a 32-year old man traveling from school in Wuhan, China to Nepal.  They used a bunch of techniques I'm not familiar with to sequence the genome.

According to the page, this sequence is published to GISAID under EPI_ISL_410301.  It seems like all of the data so far is on GISAID, so I'll treat that as my primary data source.

### [Galaxy Project](https://covid19.galaxyproject.org/)

Galaxy Project is another site centered around sharing scientific data and computational research.  On the linked page, they've put together many different workflows and examples showing how to load and analyze data.  These look at variation of the virus, how it has evolved, how different chemical components work, etc.

They also list highlights they've found from these workflows.  I didn't understand much of this, but it seems they've sussed out the useful datasets, figured out a little more how it evolves, and identified some compounds that can inhibit the spread of the virus.

From what I can tell, this site will be extremely useful.  They have lots of premade Jupyter notebooks and tools in the Anaconda cloud that I can use to load and analyze data.

## Week 2

I took at look at the COVID-19 genome with some phylogenetic trees.  Specifically, I wanted to investigate how the virus got to the US.  The notebook can be found <a href="https://nbviewer.jupyter.org/github/anderson-github-classroom/csc-448-project/blob/master/students/acengel/Project2.ipynb">here</a>.

## Week 3

I continued my efforts from last week with phylogenies.  This time, I looked at the spread of COVID-19 amongst states once it got into the US.  I had to redo some of my data scraping to stay up to date with my sources.  The notebook can be found <a href="https://nbviewer.jupyter.org/github/anderson-github-classroom/csc-448-project/blob/master/students/acengel/Project3.ipynb">here</a>.