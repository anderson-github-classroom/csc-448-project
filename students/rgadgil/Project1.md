---
jupyter:
  jupytext:
    formats: ipynb,md,py
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

## Project Part 1 (a.k.a. Project 1 on Schedule)
We don't know enough bioinformatics yet, but we can start the project by cleaning up and improving your individual websites. For Project 1, I want you to read and summarize/synthesize/discuss the following links on your own pages within the main project repo. This is definitely less of a "for grade" assignment. There is a pandemic, and I believe we can contribute meaningful bioinformatics analysis. That starts with a review of currently available resources. Please use my links as a starting point, but post more links on Slack for others. Find more links. Together we can make a contribution. Here are the links that you should investigate at a minimum.

1. I need your GitHub username, you should be able to pull, but I need to add you in order to have write permissions <a href="https://forms.gle/L3n2ytZiht11BXqU8">https://forms.gle/L3n2ytZiht11BXqU8</a>
2. https://www.gisaid.org/
3. https://foldingathome.org/
4. https://nextstrain.org/ncov
5. https://mra.asm.org/content/9/11/e00169-20
6. https://covid19.galaxyproject.org/
7. https://unanimous.ai/what-is-si/
8. https://covidtracking.com/api


### Link to clone the repository
Here is a link to the project repository.

https://github.com/anderson-github-classroom/csc-448-project

The website can be viewed at https://anderson-github-classroom.github.io/csc-448-project/.


https://www.gisaid.org/

This site contains the genetic epidimeology of COVID-19, and the genetic sequences for all strains of the virus that have been detected. This will be really useful when we learn more about gene sequences and analyzing them.


https://foldingathome.org/

A lot of diseases are caused by proteins folding incorrectly. Understanding why this happens can help us combat viruses like COVID-19 and other diseases. As Folding at Home describes it, the proteins start at an initial position and then move in various ways to create different shapes. Visualizing this requires a lot of distributed computing power, so Folding at Home requests users to let them use thie unused processing power. However, I was intrigued by the way proteins fold, and how various data points about the proteins initial and final structure can provide a glimpse into how they fold. I know a company in San Luis Obispo, called UnanimousAI, that uses something called "swarm technology" to aggregate various, seemingly insignificant data points, into larger metrics and trends. The Folding at Home website literally described folding proteins as a football game where you only know the initial positions of the players, and you have to predict the result. Unanimous AI literally does this for sports betting, where they take data points from the start of the game and use it to predict the outcome. They've had phenomal success predicting the results of such races. Here's a link to their work: https://unanimous.ai/what-is-si/


https://nextstrain.org/ncov
    
Next Strain is a very nifty data visualization dashboard that tracks the phenology of the virus (aka, it's evolution over time). You can see how the virus mutates and spreads over the course of time with various visualizations. Each strain on the dashboard provides a GISAID ID number, which allows you to view the gene sequence for that particular strain on gisaid.org


https://mra.asm.org/content/9/11/e00169-20
    
This paper is unique because it contains the complete genome sequence from someone from Wuhan, which is the origin of the virus.


https://covid19.galaxyproject.org/
    
This site contains existing open-source research and infrastructure built on the virus. The website is split intwo three sections: genomics (raw and complete sequences), evolution (positions adapted to positive vs. negative selection -- about 5 positions right now are said to be positively selected) and cheminformatics (nonstructural proteins vital to continue the life cycle of COVID)


https://covidtracking.com/api
    
This is an API endpoint I found where we can ping for COVID-19 data directly instead of having to rely on scraping, imports, etc. 

```python

```
