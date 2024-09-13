# Results

JOANA operates in two stages:

1. **Parameterization**: In the first stage, it parameterizes the significance scores from Differential Expression Analysis (DEA) using a Beta Mixture Model (BMM).
2. **Inference**: In the second stage, it estimates the probability of pathway activity within a Bayesian Network.

## Results for Step I

JOANA generates a plot demonstrating how the BMM fits the DEA results for each omics type individually.

<p align="center">
    <img src="./results/qvalues_second_moment_fitting_mixture.png" alt="BMM Fit to DEA Results" width="600"/>
</p>

Additionally, it provides a barplot showing the goodness of fit for the observed data.

<p align="center">
    <img src="./results/qvalues_second_moment_fitting_gof_hist.png" alt="Goodness of Fit" width="600"/>
</p>

## Results for Step II

### Barplots for Probabilities of Pathway Activity

JOANA offers a PDF containing barplots that display enriched pathways with probabilities exceeding 0.5 for multi-omics (cooperative) data. It also shows pathways for each single-omics modality with probabilities >= 0.5 that do not appear in the multi-omics analysis.

<p align="center">
    <img src="./results/barPlorProb.png" alt="Pathway Probabilities" width="600"/>
</p>

It also includes a graph depicting the relationships between pathways. The color of the nodes indicates the probability of pathway activity, the thickness of the edges shows the degree of interconnectedness between pathways (i.e., how many common genes they share), and the size of the nodes reflects the size of the pathways.

<p align="center">
    <img src="./results/network.png" alt="Pathway Network" width="600"/>
</p>

Furthermore, it provides box plots and beeswarm plots that illustrate the percentage of insignificant active genes (hidden-active genes) in the enriched results.

<p align="center">
    <img src="./results/percentage_beta_significant_genes_0.1_CooperativeBoxplot.png" alt="Box Plot of Insignificant Active Genes" width="600"/>
</p>

<p align="center">
    <img src="./results/percentage_beta_significant_genes_0.1_CooperativeBeesWarm.png" alt="Beeswarm Plot of Insignificant Active Genes" width="600"/>
</p>
