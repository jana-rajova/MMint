Welcome to the analysis of the Visium data for melanoma brain metastases. 

To reproduce this pipeline, SpaCET script should be run first, in order to avoid missing files as this scripts is used to estimate the celular proportions within a feature and to classify it as a tumor/stoma/interface. Following this, scripts for batch correction/transcriptome analysis and visualization scripts as well as the inferVCN can be run.

Suggested order:
1. SpaCET.Rmd
2. transcripomic_exploration.Rmd
3. decoupleR.Rmd
4. inferCNV.Rmd
5. plotting_proportions.Rmd

