<!--StartFragment-->

# **Visualization and Interpretation of Glioblastoma Gene Expression Dataset**

**Authors: Shola Paul Emiade (@Shawlyem), Abilashni Arthiswaran (@Abilashni83) and Gobikrishnan Subramaniam (@gobi3141)**

**Link to the data set:** <https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv>

**Link to the R code:** [hackbio-cancer-internship/stage 2/code.R at main · Abilashni/hackbio-cancer-internship (github.com)](https://github.com/Abilashni/hackbio-cancer-internship/blob/main/stage%202/code.R)

## **Introduction**

Glioblastoma multiforme (GBM) is the most prevalent and aggressive form of primary malignant brain tumor, particularly affecting older adults. It occurs more frequently in men than in women and is more common among Caucasians (Wu et al., 2021). GBM tumors may arise from different cell types, such as neural stem cells, astrocytes, and oligodendrocyte precursor cells, with each origin influencing the tumor’s behavior. Genetic studies have revealed several molecular subtypes of GBM, including classical, proneural, neural, and mesenchymal, which are linked to varying clinical outcomes. Furthermore, the Cancer Genome Atlas has shown that GBM tumors often harbor mutations in key pathways like p53, receptor tyrosine kinase, and the retinoblastoma pathway. While primary GBMs are often associated with mutations such as EGFR overexpression, secondary GBMs are more likely to have IDH1 and p53 mutations (Davis, 2016).



****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfhffZdg02L_upZMQuLkpC-DqdHZKwWa09ntf56ak2lwo92Fshai5Z-2sL_0dYeO7pNd9pjIV3Z4Zua6WqQ3Dg3dyt2qsOvSBxiW3OFkvNgc96q_z9H6nSOfx43_M3KfVby71P-caklh9YXcS6dmD_TfFDa?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 1: Image retrieved from Hotspot Mutations in _H3F3A_ and _IDH1_ Define Distinct Epigenetic and Biological Subgroups of Glioblastoma (**[**https://pubmed.ncbi.nlm.nih.gov/23079654/**](https://pubmed.ncbi.nlm.nih.gov/23079654/)**)**



## Results and Discussion of Heatmaps and Functional Enrichment Analysis

**Task 1:** diverging and sequential color palette for the entire dataset.

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXe5PmRxhJW8uMxjLdYUif-QrV35BksbLROTkGymWjD4_oM3JokIDpgz8yFwCQfH0C_WTFnu0xgYuV20VIGmOOEYVHJwt15ooiSRNFPwBKDjRaDwVc-pRA0r5O1ka8WrMx3UBt-P0vb4I12Qa9Uqj1pnycw?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 2: Heatmap for diverging color palette**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXf17hj330goklY3EgRNfMfR_2VzLy9-JR_dvqPHX6EzFXG2QD1NGGFBUC_3LfhQxFVVOQWlPR_KE-uG22PNW8Vn1RihnA5srCz2hzV-yk8_Y54ST93DLCfnCAcUYYON_jKjo9Vwld3IpplqjPn0psXlxzl0?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 3: Heatmap for sequential color palette**

In R programming, color palettes play an important role in visualizing data, with sequential and diverging being two common types. A sequential palette, like the default viridis in heatmaply, works best when data moves in one direction, with colors changing smoothly from light to dark to represent increasing values. On the other hand, a diverging palette is better for data with a central value, such as in a correlation matrix. It highlights both ends of the data range (e.g., negative and positive values) using different colors for the extremes and a neutral color in the middle (Galili et al., 2018).

**Task 2:** Generate variants of heatmaps 

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXekglLtkHcn7N9Fe5p-j05aUiYKHjkFu9V2P_TDapxDgRl4iLbB-kfWUTKGVBKVX-2h-9DhrNHtfZCSxfin-RR4r5kPdiW8Wlp7W3WgsMH6hG5ZoTdNvSI3CMFeCa2MWRMPQroT-XYPxdd9p3HTei_sYNk?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 4: Clustering genes (rows) alone with divergent color palette**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcqFR1OpQ7wOSRjWjXpRR-5EWoYWFQQmopfqIybyj78prXLk_zAMAugYqi39eD3DsU7vPevT-k_4gAfhT8pMUuuRFpjd1fWYUiTle0Hw4ErPsabaK0AE6_C3XO6MLCcCZjpntoxPRvjKZLadZ_O1CCrYW1w?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 5: Clustering genes (rows) alone with sequential color palette**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcLMwhKJ-VBkzhY488oz0jxqmepp9bd9b3oWZMZGRgLY91DXutozUj-W52nTazOIBD-wZUbQ1hG_Fb9CUsq1L8A0QZZebQSvg2kb2Yz-A0jDfbs4lz_wedSOn0-bHINGtDPl9o0NkeStCESrykgTtdOY3LB?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 6: Clustering samples (column) alone with divergent color palette**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXe54PEncwansKa7nfWkgj5mb_ymWi8LE3Mj5VPvA3fUFFCEwgZmJVMuPzL8yCcADQwvZVfyBHzjdOk-Y33m2PjpX6MD7iweFVOHwkgu-JcfGz0K7oTCt43AVzzZDDTjSYivvfzVl_ttCYvfQclvnEyAKABb?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 7: Clustering samples (column) alone with sequential color palette**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfHFi6WXh-HeRQt3di8bAFNsdSU2rwZW4MO0sHhSe6V4_9gYllxnBsSxh8Rw7PDlyAm1UrnKnqSN7HmDGwTTftorNgwzlShwp3yIXx9OE3Jg_ECkFTwF0ChFM7uvo39-Z2sX0ZbOemx7sO-p4s5S_QFlt5-?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 8: Clustering both genes and samples alone with divergent color palette**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcfZgxNlskNCRD1vZI_4Ibyabn6oSmzuxBhsSekg4iOGgPEn3dSrcCQfju3m8OzLwEPjVPEf9IP1o9qBOObNY4Hx8svFp7lviCQu5CF9alrvZcdLpJjAG9x9OGA52D9ZxSv47Oarn-9uJMsle9b8d_u-e8w?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 9: Clustering both genes and samples alone with sequential color palette**

Clustering rows and columns in heatmaps helps identify patterns, make complex data simpler, and show relationships. It groups similar data points, making it easier to see trends or common features, such as genes with similar expression levels. By placing related rows and columns together, clustering makes the heatmap clearer and helps interpret large datasets, making it easier to spot important patterns and connections (Engle et al., 2017).

**Task 3:** Subset genes that are significantly upregulated and downregulated

**Link for the csv file of upregulated genes:** [hackbio-cancer-internship/stage 2/data/upreg\_genes.csv at main · Abilashni/hackbio-cancer-internship (github.com)](https://github.com/Abilashni/hackbio-cancer-internship/blob/main/stage%202/data/upreg_genes.csv)

**Link for the csv file of downregulated genes:** [hackbio-cancer-internship/stage 2/data/downreg\_genes.csv at main · Abilashni/hackbio-cancer-internship (github.com)](https://github.com/Abilashni/hackbio-cancer-internship/blob/main/stage%202/data/downreg_genes.csv)

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXedgaTqwEsFlidYVnOS4w0B-knlhgigsqQHYeRvAbKX2O1usE7m5MeN0zRsHpwWDMeM55WC9UQzaW8qolfvcVa4nB729loq3epLYU8A1Rma-7YwosxsHQ3ALDN1gAmvJgyoACiSEq1hyl1dbUIm5drKkrz0?key=GW7GZ2zHVjSMNHeoiyuENQ)****

**Figure 10: Visualization of upregulated and downregulated of genes using Volcano Plot**

Transcript profiling provides important data on RNA expression changes, but detecting differentially expressed (DE) transcripts often relies on basic criteria. These include **fold change (FC)**, which shows how much expression has changed, and the **p-value**, which measures the likelihood that the change happened by chance. Using FC and p-value helps identify which transcripts show significant expression differences in an experiment (St. Laurent et al., 2013).

In this Volcano plot, the X-axis represents the **log2FoldChange** of genes. The cut off values for the  significantly downregulated genes are less than -1 and significantly upregulated genes are greater than 1. Y axis shows **-Log10 adjusted p-value** with the cut off value - 0.05. 

**Task 4:** Functional Enrichment Analysis using ShinyGO

Functional enrichment analysis (GSA) is used to identify over-represented biological annotations in a gene list and interpret molecular mechanisms and processes. The three main types of GSA are Singular Enrichment Analysis (SEA), which assesses the significance of individual annotations in a gene set; Gene Set Enrichment Analysis (GSEA), which analyzes the distribution of gene sets ranked by criteria like fold change; and Modular Enrichment Analysis (MEA), which evaluates groups of related annotations for significance by using SEA or GSEA (Garcia-Moreno et al., 2022).

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdhqxFjOAorueVfwOsx4YBlmW3TdmkZVVYacofBqAwU6WJgWdE_-_WKeFFE2VW6i3-iMxSmJHfRI-PB1kz-V9peSCXncFfR0KE-sCzBW__IjKYjJd3kwCVAGSSGQXqKj0g1I6QOkHApmHfBE65rJ9CtQDui?key=GW7GZ2zHVjSMNHeoiyuENQ)

**Figure 5: Functional Enrichment Analysis of Glioblastoma dataset**

The genes ID were submitted in ShinyGO (<http://bioinformatics.sdstate.edu/go/>)  and p-value set as 0.05. The table seems to show the results of enrichment analysis, likely from a gene set enrichment analysis (GSEA), highlighting key pathways, the number of involved genes, and various enrichment metrics. Lower FDR values indicate higher statistical significance, with pathways showing very low FDR values, such as 1.72E-08, being highly enriched. In this case, pathways like "Cytokine-cytokine receptor interaction," "Viral protein interaction with cytokine and cytokine receptor," and "Chemokine signalling pathway" have the strongest enrichment, indicated by both their very low FDR values and high fold enrichment. Hence, this table highlights immune related pathways that are highly enriched and the genes involved in these pathways play a crucial role in biological processes.

**Task 5:** Visualization of top 5 pathways

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeS41UjTK44AFsjC5AkRUUnx6XKTLARP1ef-tIxwwpvLH235WZ9hzEzqE5t2nk60N61ZGHp1fIA4A-p2gcwqGXvcx2KhHzXJ7SVga4MCjp8Ej6xGmq2ANVTlhKIIBHwei_HYHEo-w91vJ5XaUeRjHT8DhLm?key=GW7GZ2zHVjSMNHeoiyuENQ)

**Figure 6: Bar plot of top 5 enriched KEGG pathways**

**Task 6:**  Description of top 3 enriched pathways according to the biological process.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdirXoSbU9SA3_WFDbc60AGD3hKOFBO9OuE1CV0COQfCyYvrs_9v89FXx-6vTK04Gpa2nfdeZEpozaYlKsK0wPISqaST0ETqBAcHTxhPB5bk4wbXKHYgYpBtB_Ws8bWJRZsfvcc-315ihSZmzLutx3z0iVt?key=GW7GZ2zHVjSMNHeoiyuENQ)

**Figure 7: Lollipop plot according to biological process**

**Loop of Henle Development:** The loop of Henle, a U-shaped structure in the nephron, plays a critical role in reabsorbing water and sodium chloride, aiding the kidney in concentrating urine. Although the loop of Henle and glioblastoma are distinct, similar cellular pathways, such as Wnt and Notch, regulate both nephron development and glioblastoma biology (Rajakulendran et al., 2019). These pathways impact cell signaling and growth, and understanding these shared mechanisms could offer new avenues for glioblastoma treatment. Notably, there are documented cases of cerebral glioblastoma in kidney transplant recipients, a group at increased risk for cerebral tumor development due to immunosuppressive therapies (Kasiske et al., 2004).

**Proximal/ distal pattern formation:** Proximal-distal pattern formation is a crucial developmental process that organises cells and tissues along an axis from proximal (closer to the body) to distal (further from the body). In vertebrates, limb buds are often used to study this process, with genes regulating patterning through structures like the zone of polarizing activity (ZPA) and the apical ectodermal ridge (AER). Key signaling pathways, such as Hedgehog (HH), TGF-β/BMP, Wnt, and FGF, work together with genetic networks to control limb growth, positioning, and segmentation (Capdevila & Belmonte, 2001). Similarly, these pathways guide the formation of nephron segments in the kidney, from the proximal tubule to the distal tubule. In cancers like glioblastoma, disruptions in these mechanisms, particularly in pathways like Wnt and Notch, can result in uncontrolled cell growth, invasion, and improper differentiation, contributing to tumor development  (Rajakulendran et al., 2019).

**Lymphocyte chemotaxis:** Cellular interactions are facilitated by proteins called cytokines, with a large group known as chemokines responsible for mediating chemoattraction between cells (Novak et al., 2020). Lymphocyte chemotaxis refers to the movement of lymphocytes, a type of white blood cell, toward specific areas of the body in response to chemokine signals. This process is critical for immune function, allowing lymphocytes to travel to sites of infection, inflammation, or injury to trigger an immune response. Chemokine receptors on the outer cell membrane, which are part of the G-protein coupled receptor (GPCR) family, typically bind to chemokines that play a key role in cancer progression (Novak et al., 2020).



## References

Davis, M. (2016). Glioblastoma: Overview of Disease and Treatment. _Clinical Journal of Oncology Nursing_, _20_(5), S2–S8. <https://doi.org/10.1188/16.CJON.S1.2-8>.

Wu, W., Klockow, J. L., Zhang, M., Lafortune, F., Chang, E., Jin, L., Wu, Y., & Daldrup-Link, H. E. (2021). Glioblastoma multiforme (GBM): An overview of current therapies and mechanisms of resistance. In _Pharmacological Research_ (Vol. 171). Academic Press. <https://doi.org/10.1016/j.phrs.2021.105780>

Galili, T., O’Callaghan, A., Sidi, J., & Sievert, C. (2018). Heatmaply: An R package for creating interactive cluster heatmaps for online publishing. _Bioinformatics_, _34_(9), 1600–1602. <https://doi.org/10.1093/bioinformatics/btx657>

Engle, S., Whalen, S., Joshi, A., & Pollard, K. S. (2017). Unboxing cluster heatmaps. _BMC Bioinformatics_, _18_. <https://doi.org/10.1186/s12859-016-1442-6>

St. Laurent, G., Shtokalo, D., Tackett, M. R., Yang, Z., Vyatkin, Y., Milos, P. M., Seilheimer, B., McCaffrey, T. A., & Kapranov, P. (2013). On the importance of small changes in RNA expression. _Methods_, _63_(1), 18–24. <https://doi.org/10.1016/j.ymeth.2013.03.027>

Garcia-Moreno, A., López-Domínguez, R., Villatoro-García, J. A., Ramirez-Mena, A., Aparicio-Puerta, E., Hackenberg, M., Pascual-Montano, A., & Carmona-Saez, P. (2022). Functional Enrichment Analysis of Regulatory Elements. _Biomedicines_, _10_(3), 590. <https://doi.org/10.3390/biomedicines10030590>

Rajakulendran, N., Rowland, K. J., Selvadurai, H. J., Ahmadi, M., Park, N. I., Naumenko, S., Dolma, S., Ward, R. J., So, M., Lee, L., MacLeod, G., Pasiliao, C., Brandon, C., Clarke, I. D., Cusimano, M. D., Bernstein, M., Batada, N., Angers, S., & Dirks, P. B. (2019). Wnt and Notch signaling govern self-renewal and differentiation in a subset of human glioblastoma stem cells. _Genes & Development_, _33_(9–10), 498–510. <https://doi.org/10.1101/gad.321968.118>

Kasiske, B. L., Snyder, J. J., Gilbertson, D. T., & Wang, C. (2004). Cancer after Kidney Transplantation in the United States. _American Journal of Transplantation_, _4_(6), 905–913. <https://doi.org/10.1111/j.1600-6143.2004.00450.x>

Capdevila, J., & Belmonte, J. C. I. (2001). Patterning Mechanisms Controlling Vertebrate Limb Development. _Annual Review of Cell and Developmental Biology_, _17_(1), 87–132. <https://doi.org/10.1146/annurev.cellbio.17.1.87>

Novak, M., Koprivnikar Krajnc, M., Hrastar, B., Breznik, B., Majc, B., Mlinar, M., Rotter, A., Porčnik, A., Mlakar, J., Stare, K., Pestell, R. G., & Lah Turnšek, T. (2020). CCR5-Mediated Signaling is Involved in Invasion of Glioblastoma Cells in Its Microenvironment. _International Journal of Molecular Sciences_, _21_(12), 4199. <https://doi.org/10.3390/ijms21124199>



<!--EndFragment-->
