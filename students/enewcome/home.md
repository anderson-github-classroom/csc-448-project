---
permalink: /enewcome/
title: "Eric Newcomer"
excerpt: "CSC 448 Bioinformatics Algorithms Project"
---

# Project 1

**The field of bioinformatics, like other fields, relies heavily on community contributions from scientists around the world. Many of these collaborative resources can aid us in our own study of bioinformatics, and more specifically, our study of COVID-19. The ones that I have come across are listed below**.

## [GISAID Initiative](https://www.gisaid.org/)
### What is it?
The [GISAID Initiative](https://www.gisaid.org/) is a nonprofit organization that seeks to provide universal access to all influenza virus sequences, as well as their associated epidemiological and geographical data. GISAID stands for Global Initiative on Sharing All Influenza Data, and is a partnership between the German government and a registered nonprofit called Freunde of Gisaid (Friends of Gisaid). In addition to providing an open-access, influenza virus database, they provide training workshops around the world and aid in the development of new tools for GISAID data analysis.

### Why is it important?
By effectively "open-sourcing" influenza virus sequence data, GISAID has been able to lift the hurdles and restrictions that would otherwise stop bioinformaticians, scientists, and students around the world from finding useful insights about the virus data that they can then provide to the scientific community. It is an extremely important resource to scientists and bioinformaticians alike, especially as we try to learn more about COVID-19.  

### Other cool things/facts:
- In order to access the free database, you must first submit an application that includes basic personal information and an agreement to their database terms. After an applicant's identity is confirmed, they will receive access credentials to the database.

## [Folding@Home](https://foldingathome.org/)
### What is it?
Proteins are essential for keeping humans healthy and are assembled when they fold. However, protein misfolding has serious consequences that can lead to the development of many diseases. Folding@Home is a community-driven, distributed computing project that aims to analyze mass amounts of protein folding simulation data. The simulations are time-intensive, which is why the project relies on the contributions of citizen-scientists around the world, who volunteer to run protein dynamics simulations on their personal computers. Users contribute by downloading the Folding@Home software, which runs simulations while you are doing other things.

### Why is it important?
Folding@Home is an amazing project that allows anyone with a computer and an internet connection to contribute to disease research. There are a massive amount of computations that need to be done to study protein dynamics; by sharing unused computer power, more insights can be gained and more cures can be researched. 

### Other cool things/facts:
- The Folding@Home project has a [Diseases page](https://foldingathome.org/diseases/) that includes more information as well as a timeline of their research findings *for every disease that they have studied*. 
- All of the diseases they study are organized into 3 categories, being cancer, infectious diseases, and neurological diseases. Their research encompasses a wide range of diseases, such as the Zika virus, Ebola virus, Alzheimer's disease, and Parkinson's disease.

## [Next Strain](https://nextstrain.org/ncov/global)
### What is it?
Next Strain is an open-source project that seeks to combine, understand, visualize, and track pathogen genome data in real time. The primary goal of the project is to increase epidemiological understanding and improve the community's response to the outbreak of viral diseases. Their latest [data and analysis page](https://nextstrain.org/ncov/global) contains a dashboard that analyzes and displays the Genetic epidemiology of the novel coronavirus. The main visualization shows the phylogeny (evolutionary history) of the virus, which can be displayed globally and by country.

### Why is it important?
Data visualization in general is an extremely powerful tool that allows for increased comprehension of data analysis, since you don't (usually) need to have scientific/statistical background to understand what is happening in a picture. Furthermore, by visualizing the data in a multitude of ways (similar to the different tree options on the NextStrain dashboard) more relationships/correlations/insights can be found at the click of a button. The Next Strain data visualization dashboard allows users to view how pathogens change over time in a plethora of ways, enabling people everywhere to see and try to understand how viruses  evolve.

### Other cool things/facts:
- This project is made possible using data from GISAID!
- The dashboard provides analysis for 10 different pathogens, including ebola, ncov (novel coronavirus), and zika.
- Next Strain provides an open-source toolkit that allows anyone to create visualizations like the ones that are seen on the dashboard.

## [Complete Genome Sequence of a 2019 Novel Coronavirus (SARS-CoV-2) Strain Isolated in Nepal](https://mra.asm.org/content/9/11/e00169-20)
### What is it?
The resource listed above is a link to an article by the American Society for Microbiology that contains the complete genome sequence of a COVID-19 strain isolated in Nepal. The sequence was obtained via an oropharyngeal swab specimen of a 32 year-old Nepalese student with coronavirus disease 2019 (COVID-19), who had returned to Nepal after traveling to Wuhan, China. The student tested positive by real-time reverse transcriptase PCR developed by the University of Hong Kong. Full genome comparison revealed that there was a >99.99% match with two previously sequenced genomes from GenBank. The SARS-CoV-2 genome contains a single, positive-stranded RNA that is 29,811 nucleotides long, being 29.86% adenosines, 18.39% cytosines, 19.63% guanines, and 32.12% thymines.

### Why is it important?
The article and its findings are extremely important to the bioinformatics and scientific community. It provides a complete genome sequence of the COVID-19 virus, which verifies previous findings and allows for scientists and bioinformaticians alike to start analyzing the genome.

### Other cool things/facts:
- The patient was a 32 year-old Nepalese man who attended the Wuhan University of Technology as a student.  
- His symptoms were cough, mild fever, and throat congestion.

## [Best practices for the analysis of SARS-CoV-2 data: Genomics, Evolution, and Cheminformatics](https://covid19.galaxyproject.org/)
### What is it?
The above website is a project created by the [Galaxy Project](https://galaxyproject.org/) that seeks to provide the tools and infrastructure needed for users to analyze viral datasets. The three analysis types that are featured are Genomics, Evolution, and Cheminformatics. What each analysis type encompasses is listed below:

**Genomics**:
- Assembly
- MRCA timing
- Variation analysis
- Selection adn recombination

**Evolution**:
- Natural Selection Analysis
- Visualizations
- Observable Notebooks

**Cheminformatics**:
- Compound enumeration
- Generation of 3D conformations
- Docking
- Scoring
- Selection of compounds for synthesis

### Why is it important?
This project, as well as the other tools provided by Galaxy, allows everybody to access powerful data science tools that enables them to participate and learn about science/data analysis. The tools are not only for computational research; they also serve as free, science/data analysis education for the public.

### Other cool things/facts:
- The Genomics section relies on data provided by GISAID!
- There are 397 sites showing intra-host variation across 33 samples (with frequencies between 5% and 95%). (Genomics section)
- At present, ~5 genomic positions may merit further investigation because they may be subject to diversifying positive selection. (Evolution section)
- 40,000 compounds were considered to be likely to bind, and were chosen based on recently published X-ray crystal structures; 500 high scoring compounds were identified. (Cheminformatics section)

---

# Project 2
Whenever there is a outbreak due to a new virus, scientists and researchers must determine where the virus came from, when it appeared, and how it was transmitted. To do this, they use evolutionary trees, or phylogenies. For this week's project, I wanted to construct phylogenies to answer the following questions: **What is the variation in the most uncommon sequences in the most affected countries?** The data analysis and phylogenetic trees can be found [here](../Project2.ipynb).








```python

```
