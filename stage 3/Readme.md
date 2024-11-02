
**Biomarker Discovery through Machine Learning and Differential Expression Analysis of Stomach Adenocarcinoma**

**Authors**: Olabode Kaosara (@TheResearchNerd)

Ala Khaled (@Alaa043)

Ahmed Abdel-Masqoud (@Abdel-Maqsoud)

Obianujunwa Martha (@obianujunwa)

Raqeebah Rafiu (@Raqeebahh)

Niguse Kekile Lema (@King)

Abilashni Arthiswaran (@Abilashni83)

Benedict Orile Ajunku (@orile)

**INTRODUCTION**

Stomach adenocarcinoma, or gastric cancer, is a formidable malignancy that poses a significant threat to the digestive system worldwide. \[1]. It ranks fourth in incidence after lung, breast, and colorectal cancer and second in cancer-related mortality after lung cancer. One million people are diagnosed with gastric cancer worldwide, primarily in developing countries \[2]

![STAD](https://github.com/user-attachments/assets/00d79f3b-f238-48a3-b716-c15d6af7b21b)
\[3]

**BRIEF DESCRIPTION OF DATA**

To conduct our analysis, we employed the STAD dataset from TCGA, which contains the clinical and omics data. Using the TCGAbiolinks package in Bioconductor, we leveraged the GDCdownload and GDCprepare functions to retrieve RNA-seq data in RStudio. We subsequently streamlined the metadata dataset by selecting relevant features and reduced the sample pool to an equal representation of 20 female and 20 male participants. To ensure data accuracy, we normalized the raw expression counts based on gene length and eliminated genes exhibiting low expression levels.

**METHODOLOGY**

To identify differentially expressed genes (DEGs) between the two groups, we conducted a differential expression analysis (DEA) using the TCGAbiolinks library and the edgeR package. Genes were considered differentially expressed if they exhibited an absolute log fold change greater than 1 and an adjusted p-value less than 0.05. This analysis yielded 654 upregulated and 926 downregulated genes.

To gain insights into the biological processes associated with these DEGs, we performed an enrichment analysis using the TCGAbiolinks library. The K-nearest neighbors (KNN) model, a machine learning algorithm commonly used for classifying biological data, was employed to predict gender based on gene expression profiles.

**RESULT AND DISCUSSION**

The results of this analysis were visualized in a bar plot. From the output of the KNN classification, it can be concluded that the model was able to predict the validation set with an accuracy of 92%.

This study successfully integrated machine learning and differential expression analysis to explore potential biomarkers in stomach adenocarcinoma. While differential expression analysis identified key upregulated and downregulated genes, the K-Nearest Neighbor model's suboptimal performance on the test dataset suggests that larger datasets and model refinement are needed to enhance predictive accuracy. It is recommended that further research should focus on expanding the dataset and exploring more sophisticated machine-learning algorithms to improve biomarker identification and ultimately aid in the early detection and treatment of stomach adenocarcinoma.
<img width="499" alt="HEATMAP F VS M WITH DENDOGRAM AND LEGEND" src="https://github.com/user-attachments/assets/5ab686f1-c90d-4959-83c5-df453e9bad58">

<img width="499" alt="HEATMAP F VS M WITHOUT DENDOGRAMS" src="https://github.com/user-attachments/assets/44f9da4d-4fee-4db0-a947-357196644a9c">

![dnreg genes](https://github.com/user-attachments/assets/72722373-a225-4dfc-9c51-84587f494e45)

![upreg genes](https://github.com/user-attachments/assets/50905db9-af82-43a3-ae52-1f8a3c0210c6)

<img width="499" alt="VOLCANO F VS M" src="https://github.com/user-attachments/assets/1497c922-f2df-4fe3-94a2-efcc1e0b3c6b">

**BOXPLOT**

<img width="506" alt="BOXPLOT" src="https://github.com/user-attachments/assets/277b2096-ceb1-4866-a495-fe46465d9378">

**REFERENCES**

<!--[if !supportLists]-->1.     <!--[endif]-->Hermann, P.C. _et al._ (2007). Distinct populations of cancer stem cells determine tumor growth and metastatic activity in human pancreatic cancer. _Cell Stem Cell._

<!--[if !supportLists]-->2.     <!--[endif]-->Miraglia, S. _et al._ (1997). A novel five-transmembrane hematopoietic stem cell antigen: isolation, characterization, and molecular cloning. _Blood_.

<!--[if !supportLists]-->3.     <!--[endif]-->Ji Min Choi, Stanley Oiseth, Lindsay Jones, Evelin Maza. (2022). \_Gastric cancer\_. https\://www\.lecturio.com/concepts/gastric-cancer/
