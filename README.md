# GWAS Analysis of Lingering Ash
The goal of this analysis is to identify genomic regions that may harbor some genetic mechanism to EAB resistance in ash (*Fraxinus*). 

<img width="500" alt="Genomics guides" src="https://github.com/user-attachments/assets/c20df508-87c1-42b7-ac27-547db78bc1c2" />

## Overview of Input Data
I have 369 trees ranging from pure green ash (*F. pennsylvanica*) to pure white as (*F. americana*) from the USFS NRS with both genomic and phenotypic data. There is also a gradient of hybrids. These trees have been assessed for disease phenotype which is a proxy for resistance.
### Distribution of Families
<img height="500" alt="histograms of the number of families, mothers, and fathers" src="https://github.com/user-attachments/assets/1d5366c8-fcc2-481d-93b2-32b8bc929d12" />
Kinship will need to be accounted for in the GWAS models. I have individuals from 24 families. These are lingering x lingering crosses.


## Step 1: Phenotype Analysis
I was sent the phenotype information from the collaborators at the USFS. There are a few values that I think will be useful to determine phenotype: pHK and pL34. pHK is the proportion of larvae that were killed by the host tree. pL34 are the proportion of larvae that were at instars 3 or 4. In the heritability paper, it mentions that pHK is the best metric for phenotype for disease resistance/susceptibility. I ran a few basic linear models to determine association between the phenotype variables (below). This was ran in R and the associated script is in folder **Step1**.

<details>
  <summary>Phenotype Visuals</summary>

<img width="850" height="768" alt="linear_models" src="https://github.com/user-attachments/assets/cf24bbe4-46e5-4df6-ab8e-6995e886a4d9" />

As evidenced by the linear models, pHK and pL34 have the most significant relationships. Then, I checked to see how these values are distributed across families.

<img height="500" alt="mean and 95% CI plots of pHK and pL34 per family" src="https://github.com/user-attachments/assets/0a676705-e17d-4490-af9e-e3c1e284a580" />

Some families clearly score better than others. We want a high pHK and a low pL34. This means the tree is killing off the larvae before they reach instars 3 and 4.
</details>

# Attempt 1: 369 samples (both pheno and geno data)
## Step 2: Variant Calling
For this, I am using the **vary_cool** pipeline developed by the incredible, amazing, wonderful Staton lab. The github is available [here](https://github.com/statonlab/vary_cool). All scripts used in this are in folder **Step2**. It takes raw fastq files and gives you a nice VCF file. 

  1) I moved all of the files from the shared directory to my scratch folder with the script *00.move_files.sh*.
  2) I used Claude to write a script that takes the sampleID from a csv file and moves it into a folder. In this case, I made a csv file with all of the individuals that we have phenotype data for and directed them into a new folder. This is script *00.move_files_csv.sh*. This does not use slurm, so execute it [bash organize_fastq.sh ./fastq_files samples.csv ./matched_samples].
  3) Next, I made a metadata table. Vary_cool is vary_particular about the metadata file so I copied some script from Zane. It is in *01.create.metadata*. That metadata file needs to be in the parent GWAS folder so it can be accessed by nextflow.
  4) Finally, I ran the script. It took forever to make it to the top of the queue. The script is *02.run_vary_cool.sh*.
