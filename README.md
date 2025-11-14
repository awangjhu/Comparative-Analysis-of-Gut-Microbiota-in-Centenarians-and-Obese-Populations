# Check-in 2: Comparative Analysis of Gut Microbiota in Centenarians and Young Populations
**Alan Wang and Carsten Hoeke**

---

## Addressing Prior Feedback
Since the last check-in, our focus has been building a functional and reproducible metagenomic analysis pipeline to handle large data sets. However, as mentioned in the previous check-in, we have had issues with our workflow related to Kraken2 and Braken database compatibility (Braken processing), which prevented us from generating accurate abundance estimates. To address this, we have been working to identify what is causing these issues and redesigned our directory structure to create a Kraken2 to Braken pipeline. In addition, our previous check-in 1 readme.md did not contain our scripts to generate the figure as well as the figures inside of the readme.md file ended up crashing. Since then, we have addressed both issues and check-in 1 readme.md is fixed.

---

## Progress
After implementing the previous feedback, we began streamlining our workflow by improving the clarity and usability of our data-processing scripts. We added comments to each major code block to explain its role in the overall pipeline, making the script easier for new users to navigate. We also replaced long, hard-coded file paths with clean variable definitions, which improves readability and ensures the script can run smoothly across different systems and directory organizations. Next, we continued the development of our Braken statistical analysis using a test batch with a singular data set from the centenarian and young cohorts. The original paper includes a third cohort of “elderly” which was also included in the Braken analysis.1 This novel cohort was initially processed using the Kraken pipeline, and subsequently ran through the Braken analysis. The same abundance composition at the Class taxonomic level was conducted between the three cohorts, and yielded results consistent with the less rigorous Kraken analysis. To further elucidate potential factors with functional relevance to the differential composition, we first needed to detect differentially abundant taxa. The differential abundance analysis was performed using the ANCOM-BC2 R package, which contains a series of functions that output a list of enriched, class-level taxa. This analysis was restricted to the centenarian and young cohorts using our preliminary single dataset for pipeline construction. These data were plotted in a volcano-plot, where logFC > 0 represents higher enrichment in the young cohort and logFC < 0 represents higher enrichment in the centenarian cohort. The top hits were Mollicutes (Higher in Centenarians, logFC = -0.822), Coriobacteriia (Higher in Young, logFC = 0.675), Bacilli (Higher in Centenarians, lgFC = -0.0.641).

![Figure 1](figures/Braken%20Stacked%20Bar%20Figure.png)  
**Figure 1.** Bacterial composition from Braken pipeline.

![Figure 2](figures/brckn_ancom_volcano.png)  
**Figure 2.** Diffenetial abundance analysis of bacterial classes between young and centenarian samples.

Although our current dataset is still limited to one sample per group, our findings show age-associated trends that seem quite consistent with findings from the published observations. Notably, there seems to be a decline in Erysipelotrichia, Negativicutes, and Alphaproteobacteria, while an increase in Verrucomicrobiae, Bacilli, Bacteroidia (Figure 1). We also used ANCOM-BC2 on our Braken class-level abundance data to determine differential taxa between young and centenarian samples. Some classes showed fold changes; however, no taxa reached statistical significance because we used a sample size of 1 per group (Figure 2). 

In the future, we plan to expand the cohort to enable formal statistical testing with greater power, but our current results confirm that the differential abundance pipeline is functioning as intended. We would also like to include the elderly cohort and create another method to visualize the changing abundance over time (pseudo-time?). Beyond improving sample size, our next goal is to connect taxonomic differences to potential functional implications. To accomplish this, we will use the BugSigDB R package, which links taxa to microbiome signatures from the literature. This approach will help identify whether the classes that differ between centenarians and young individuals are associated with specific biological pathways or host conditions, allowing us to infer functional distinctions in the gut microbiome between the two age groups.

---

## Project Organization

![Updated BUDDY](Screenshot%202025-11-13%20225732.png)

Our BUDDY repository now uses Kraken2, Braken, and downstream merging steps. All our Braken outputs are stored under results/braken/, and merged relative-abundance tables for each taxonomic level are found in results/merged/. We can then use these outputs into R for visualization and other types of statistical analysis/visualization.

---

## Challenges and Questions
One key limitation in the BUDDY project is the extremely small sample size. With only one sample per group, no statistical testing is possible, and observed differences can only be characterized as qualitative. Fortunately, we have validated multiple analysis pipelines that are functional and robust. As a result, it is simple and possible to run multiple samples using our pipeline to generate reliable data. Increasing the number of samples to around 10-20 per group remains, arguably, the most important step in the BUDDY project. In addition, we plan on building an R workflow for fold-change visualization and taxonomic trajectories to better understand what specific class of bacteria are changing due to age, such as in other publications that look at gut microbiota.2 

A couple of outstanding questions remain. Would it be possible to implement these scripts in Rockfish to make scaling this study more feasible? Currently, with between two and seven datasets, we are working with over 150G of data. To expand this to all of the available data with our current machinery is nigh impossible. Secondly, if we wanted to include the elderly cohort in our ANCOMB differential abundance pipeline, would we create three separate volcano plots (E vs. Y, Y vs. C, C vs. E)?

---

## References
1 Wu, L. et al. A Cross-Sectional Study of Compositional and Functional Profiles of Gut Microbiota in Sardinian Centenarians. mSystems 4, 10.1128/msystems.00325-00319 (2019). https://doi.org:10.1128/msystems.00325-19  
2 Rosell-Díaz, M. et al. Gut microbiota links to serum ferritin and cognition. Gut Microbes 15, 2290318 (2023). https://doi.org:10.1080/19490976.2023.2290318
