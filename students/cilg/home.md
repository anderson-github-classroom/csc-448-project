---
permalink: /cilg/
title: "Bart Ilg"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---
Week 1
------

GISAID
------
GISAID is short for Global Initiative on Sharing All Influenza Data, and involves both public and private partnerships between a non-profit organization, the German government, The Singaporean government, and the United States government. GISAID is hosted by Germany, and is supported by many philanthropists.


Folding at Home
---------------
Folding at Home is an organization who provides software that allows users to contribute their own computing power to help solve complex disease research problems. Users can form groups, and the project has been endorsed by many large companies like intel, nVIDIA, and AMD. Many other influencers have endorsed the project and created their own groups, including Linus Tech Tips, a popular tech centric youtube channel with over 10 million subscribers, and Mark Kern, the team lead for the original release of the popular video game World of Warcraft.


Nextstrain
----------
Nextstrain is a service that provides large amounts of data relating to viruses. In this case, their information on the COVID-19 virus is of interest. They have data collected from around the globe, and have graphs and maps to help track the spread of the disease.


American Society for Microbiology
---------------------------------
The American Society for Microbiologyreports that a complete genome sequence for the CVOID-19 virus has been found in Nepal. This genome could prove very useful, and it supplies researchers with information on its composition. using this genome, we can use tools provided by other services, such as Folding at Home, in order to learn more about the disease and help researchers try to develop a vaccine.


Galaxy Project
--------------
Galaxy Project provides many resources to help research COVID-19, ranging from data on its genome, to evolutionary data. With this information researchers can try to determine how it behaves, and when it has evolved to adapt to new environments.


Week 2
------
When I began working on part 2 of the project, I was unsure what to do, so I decided to experiment with the software virulign. Virulign, can be found at https://github.com/rega-cev/virulign . I sought to test virulign by comparing their provided SARS-CoV-2 Genome to another copy of the genome. SARS-CoV-2 is the name of the virus itself, while COVID-19 is the name of the disease it causes. The other copy I decided to compare it to was provided by Galaxy Project. Galaxy project's COVID-19 Genome can be found here:  https://covid19.galaxyproject.org/genomics/4-Variation/current_complete_ncov_genomes.fasta . The interesting thing that I found, was that one part of the sequence had an error. I'm interested to see what this could mean, and intend to look further into it during week 3. Could one of the sequences be wrong? If thats the case, I wonder if there are researchers out there using incorrect data. COuld this be a difference between 2 strains? I have heard that there is supposed to be a stronger and weaker strain of the virus. Could it be an issue with virulign? If so, are there more accurate alternatives out there? I intend to find out why this happened, and further expand upon the issue during week 3.


Week 3
------
Upon further inspection of virulign data, as well as a bit of research, I cam eot an understanding of what the output meant. Having little experience with this type of data, an error arising in the comparison of what should be the same 2 genomes, stood out as a potential problem with either one of the genomes, or a bug in the software. However, further research clarified that this isn't necessarily the case. In the virulign data , some of the outputs that were successful also came paired with a combination of letters and numbers. during my research I discovered that these letter number combinations in the comparison represented specific mutations. These mutations mean that while both of the presented genomes are of the same virus, they aren't completely identical. With enough mutation they could even be considered separate strains. This is likely where the alignment error came from as well. A mutation caused a slight different between the genomes, but they are both still SARS-COV-2. This information is relieving, showing that researchers aren't being given conflicting data, and virulign should be working properly.