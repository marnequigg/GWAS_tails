# GWAS Analysis of Lingering Ash
The goal of this analysis is to identify genomic regions that may harbor some genetic mechanism to EAB resistance in ash (*Fraxinus*). 

<img width="500" alt="Graphical overview of GWAS project" src="https://github.com/user-attachments/assets/5e5ebc92-f63a-4be1-9cca-414a49563b4f" />

## Overview of Input Data
I have 369 trees ranging from pure green ash (*F. pennsylvanica*) to pure white as (*F. americana*) from the USFS NRS with both genomic and phenotypic data. There is also a gradient of hybrids. These trees have been assessed for disease phenotype which is a proxy for resistance.
### Distribution of Families
<img height="500" alt="histograms of the number of families, mothers, and fathers" src="https://github.com/user-attachments/assets/1d5366c8-fcc2-481d-93b2-32b8bc929d12" />
Kinship will need to be accounted for in the GWAS models. I have individuals from 24 families. These are lingering x lingering crosses.


## Step 1: Phenotype Analysis
I was sent the phenotype information from the collaborators at the USFS. There are a few values that I think will be useful to determine phenotype: pHK and pL34. pHK is the proportion of larvae that were killed by the host tree. pL34 are the proportion of larvae that were at instars 3 or 4. In the heritability paper, it mentions that pHK is the best metric for phenotype for disease resistance/susceptibility. I ran a few basic linear models to determine association between the phenotype variables (below).

<img width="850" height="768" alt="linear_models" src="https://github.com/user-attachments/assets/cf24bbe4-46e5-4df6-ab8e-6995e886a4d9" />

As evidenced by the linear models, pHK and pL34 have the most significant relationships. Then, I checked to see how these values are distributed across families.

<img height="500" alt="mean and 95% CI plots of pHK and pL34 per family" src="https://github.com/user-attachments/assets/0a676705-e17d-4490-af9e-e3c1e284a580" />

Some families clearly score better than others.
