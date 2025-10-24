# Comparative Analysis of Gut Microbiota in Centenarians and Obese Populations

## Description
In the landmark paper “Healthspan and lifespan extension by fecal
microbiota transplantation into progeroid mice” (Bárcena et al., Nature Medicine, 2019), specific microbiota composition trends, such as an increased abundance of pro-inflammatory Proteobacteria and Cyanobacteria, are linked to gut dysbiosis and accelerated aging. Similarly, it was discovered that centenarians also saw an abundance of specific microbial species like Akkermansia muciniphila. In contrast to the centenarians, other studies have also concluded that those with lower healthspan, such as obese American adults, are significantly enriched with pro-inflammatory communities while containing lower butyrate-producing communities (Gordon et al., Nature, 2009). As a result, there is strong evidence to believe certain microbial communities are linked to both lifespan and healthspan. In this study, we will use 16S rRNA analysis to identify microbiota taxa that determine longevity between Sardinian centenarians and obese American populations.

## Example Image
![Example Image](https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41591-019-0504-5/MediaObjects/41591_2019_504_Fig1_HTML.png?as=webp)

## Data
- Sardinian Centaurian: https://www.ebi.ac.uk/ena/browser/view/PRJEB25514  
- Obese American: https://www.ebi.ac.uk/ena/browser/view/PRJEB11419  
  *(note: although this contains BMI in metadata, we will parse out participants that have a BMI over 30)*

## Software
- Kracken (for raw data): https://ccb.jhu.edu/software/kraken/
- ANCOMBC pipeline (for differential analysis): https://github.com/FrederickHuangLin/ANCOMBC  
- BugSigDB package (for GO-like analysis): https://bugsigdb.org/Main_Page  
- SPARCC package (for co-occurrence network): https://github.com/dlegor/SparCC  

## Proposed Steps

### < 10-Hour Goal
We will begin by conducting a quantitative taxonomic analysis at the phylum level of the 16S rRNA data sets, creating a stacked barchart of the composition of the gut flora between centenarian and obese cohorts. Statistical tests will be applied to determine whether statistical differences exist.

### Stretch Goal 1
Perform a “GO” like enrichment analysis of the taxa to draw conclusions on abundant or sparse microbial populations within each data set and their association with drug resistance, diets, age, and disease.

### Stretch Goal 2
Perform a co-occurrence networking analysis using a pre-built pipeline to understand whether microbial communities within a cohort are connected or fragmented.
