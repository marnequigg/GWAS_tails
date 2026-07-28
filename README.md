# GWAS Analysis of Lingering Ash
The goal of this analysis is to identify genomic regions that may harbor some genetic mechanism to EAB resistance in ash (*Fraxinus*). 
<details>
  <summary>Graphical abstract for study</summary>
<img width="850" alt="Genomics guides" src="https://github.com/user-attachments/assets/c20df508-87c1-42b7-ac27-547db78bc1c2" />
</details>

## Overview of Input Data
I have 369 trees ranging from pure green ash (*F. pennsylvanica*) to pure white as (*F. americana*) from the USFS NRS with both genomic and phenotypic data. There is also a gradient of hybrids. These trees have been assessed for disease phenotype which is a proxy for resistance.
### Distribution of Families
<img height="500" alt="histograms of the number of families, mothers, and fathers" src="https://github.com/user-attachments/assets/1d5366c8-fcc2-481d-93b2-32b8bc929d12" />
Kinship will need to be accounted for in the GWAS models. I have individuals from 24 families. These are lingering x lingering crosses.


## Step 1: Phenotype Analysis
I was sent the phenotype information from the collaborators at the USFS. There are a few values that I think will be useful to determine phenotype: pHK and pL34. pHK is the proportion of larvae that were killed by the host tree. pL34 are the proportion of larvae that were at instars 3 or 4. In the heritability paper, it mentions that pHK is the best metric for phenotype for disease resistance/susceptibility. I ran a few basic linear models to determine association between the phenotype variables (below). This was ran in R and the associated script is in folder **Step1**.

<details>
  <summary>Script</summary>
---
title: "phenotype_prelim"
author: "m quigg"
date: "2026-07-14"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Phenotype Preliminary Analysis
The purpose of this is to take a look at the resistance phenotype data the the USFS provided. They inoculated the trees with EAB eggs, then looked at the number of larvae that were remaining. 

## Goal 
See the distribution of the data to get an idea of what would be good covariables for the GWAS.

## Prep
```{r, include=FALSE}
###activate packages
library(tidyverse)
library(patchwork)
library(AICcmodavg)

###set working directory
setwd("C:/Users/mquig/Desktop/ash/GWAS_tails/phenotype")

###import the file
pheno_data=read.csv("LA-ProgenyTallyQCSAS_71426.csv", header = T)
```

After chatting with Zane, we realized that I was looking at all individuals, not just the ones we sequenced. So now, I want to rerun everything and incorporate it.
```{r}
###import the file
inds_seq=read.csv("USFS_Tails_Metadata_2023.xlsx - Sex_IDs_11.22.25.csv", header = T)

###redo the row names
inds_seq=inds_seq %>% 
  rename("IndShortName"=X.rname)

###join the dataframes
pheno_data=left_join(inds_seq, pheno_data, by = "IndShortName")

###remove the inds that we have sequence data for but no phenotype data
pheno_data=pheno_data %>% 
  drop_na(pL34)
```


## USFS Relevant Stats

### Proportion Host Kill (pHK or pTK)
This is the number of larvae that host killed divided by the number of good eggs.
```{r}
###basic stats
##max killed?
max(pheno_data$pHK)
#0.9166667
##min killed
min(pheno_data$pHK)
#0=none

##average
mean(pheno_data$pHK)
#0.2391576

###distribution
hist(pheno_data$pHK)
#right skewed

pHK_plot=ggplot(pheno_data, aes(x = pHK)) +
  geom_histogram(bins = 10, color = "magenta4", fill = "pink") +
  labs(title = "Proportion of Host Kills") +
  theme_classic()
pHK_plot
```
Ok so hist() broke it into bins with 0.1 intervals. There are a lot of individuals in the 0-0.1 pHK category.

### Realized Dose
This is the number of good eggs per tree.
```{r}
max(pheno_data$RlzDose)
min(pheno_data$RlzDose)
hist(pheno_data$RlzDose)
RlzDose_plot=ggplot(pheno_data, aes(x = RlzDose)) +
  geom_histogram(bins = 10, color = "aquamarine3", fill = "aquamarine") +
  labs(title = "Realized Dose of Good Eggs") +
  theme_classic()
RlzDose_plot
```

### Proportion 3 and 4 Instar Larvae
This is the proportion of recovered larvae that made it to instars 3 and 4. This is a proxy for EAB survival.
```{r}
max(pheno_data$pL34)
min(pheno_data$pL34)
hist(pheno_data$pL34)
pL34_plot=ggplot(pheno_data, aes(x = pL34)) +
  geom_histogram(bins = 10, color = "purple4", fill = "mediumpurple1") +
  labs(title = "Proportion of Larvae at Instars 3 or 4") +
  theme_classic()
pL34_plot
```

### Large Larvae
This is the number of large, surviving larvae at instars 3 and 4.
```{r}
max(pheno_data$LgLarv)
min(pheno_data$LgLarv)
hist(pheno_data$LgLarv)
LgLarv_plot=ggplot(pheno_data, aes(x = LgLarv)) +
  geom_histogram(bins = 10, color = "yellowgreen", fill = "darkolivegreen1") +
  labs(title = "Number of Large Larvae") +
  theme_classic()
LgLarv_plot
```
### Make a figure
```{r}
(pHK_plot + pL34_plot)/(RlzDose_plot + LgLarv_plot)
```
### Comparisons
Meg suggested that I run some little analyses to see if any of these variables are associated.

#### pHK and pL34
```{r}
###pHK vs pL34
##make the graph
pHK_vs_pL34_plot=ggplot(pheno_data, aes(x = pHK, y = pL34)) +
  geom_point() +
  geom_smooth(method = "lm", color = "pink", fill = "magenta4") +
  theme_classic()
pHK_vs_pL34_plot

##make the model
pHK_vs_pL34_lm=lm(pL34 ~ pHK, pheno_data)
#intercept=0.8062
#-0.8184


###pL34 vs pHK
##make the graph
pL34_vs_pHK_plot=ggplot(pheno_data, aes(x = pL34, y = pHK)) +
  geom_point() +
  geom_smooth(method = "lm", color = "aquamarine3", fill = "aquamarine") +
  theme_classic()
pL34_vs_pHK_plot

##make the model
pL34_vs_pHK_lm=lm(pHK ~ pL34, pheno_data)
#intercept=0.7851
#-0.8748
```

### pHK and RlzDose
```{r}
###pHK vs RlzDose
##make the graph
pHK_vs_RlzDose_plot=ggplot(pheno_data, aes(x = pHK, y = RlzDose)) +
  geom_point() +
  geom_smooth(method = "lm", color = "pink", fill = "magenta4") +
  theme_classic()
pHK_vs_RlzDose_plot

##make the model
pHK_vs_RlzDose_lm=lm(RlzDose ~ pHK, pheno_data)
pHK_vs_RlzDose_lm
#intercept=388.11
#33.58

###RlzDose vs pHK
##make the graph
RlzDose_vs_pHK_plot=ggplot(pheno_data, aes(x = RlzDose, y = pHK)) +
  geom_point() +
  geom_smooth(method = "lm", color = "purple4", fill = "mediumpurple1") +
  theme_classic()
RlzDose_vs_pHK_plot

##make the model
Rlzdose_vs_pHK_lm=lm(pHK ~ RlzDose, pheno_data)
Rlzdose_vs_pHK_lm
#intercept=-0.137213
#0.001053
```

#### RlzDose and pL34
```{r}
###pHK vs pL34
##make the graph
RlzDose_vs_pL34_plot=ggplot(pheno_data, aes(x = RlzDose, y = pL34)) +
  geom_point() +
  geom_smooth(method = "lm", color = "purple4", fill = "mediumpurple1") +
  theme_classic()
RlzDose_vs_pL34_plot

##make the model
RlzDose_vs_pL34_lm=lm(pL34 ~ RlzDose, pheno_data)
RlzDose_vs_pL34_lm
#intercept=1.000288
#-0.001067


###pL34 vs pHK
##make the graph
pL34_vs_RlzDose_plot=ggplot(pheno_data, aes(x = pL34, y = RlzDose)) +
  geom_point() +
  geom_smooth(method = "lm", color = "aquamarine3", fill = "aquamarine") +
  theme_classic()
pL34_vs_RlzDose_plot

##make the model
pL34_vs_RlzDose_lm=lm(RlzDose ~ pL34, pheno_data)
pL34_vs_RlzDose_lm
#intercept=418.52
#-36.39
```
### compare AIC scores
```{r}
AIC(pHK_vs_pL34_lm, pHK_vs_RlzDose_lm, pL34_vs_pHK_lm, pL34_vs_RlzDose_lm, Rlzdose_vs_pHK_lm, RlzDose_vs_pL34_lm)
```
I don't think this is really telling us anything.

### Make a figure
```{r}
(pHK_vs_pL34_plot + pHK_vs_RlzDose_plot) / (pL34_vs_pHK_plot + pL34_vs_RlzDose_plot) / (RlzDose_vs_pHK_plot + RlzDose_vs_pL34_plot)
```



## Other Relevant Visuals
I want to make a stacked bar chart with the good and bad eggs per tree. This is going to have a lot of bars, but I think it'll work.
```{r}
###create the data frame
egg=pheno_data %>% 
  dplyr::select(IndShortName, GoodEgg, BADEGG) %>% 
  rename("Good" = GoodEgg, "Bad" = BADEGG)
egg_long=pivot_longer(egg, cols = c(Good, Bad), names_to = "egg", values_to = "count")

###graph it
ggplot(egg_long, aes(x = IndShortName, y = count, fill = egg)) +
  geom_bar(position="stack", stat="identity")
```
So that doesn't look like what I wanted, but it works!


### Parents
```{r}
summary_families=pheno_data %>% 
  group_by(Family) %>% 
  summarise(n=length(Family)) %>% 
  ungroup()
summary_moms=pheno_data %>% 
  group_by(FemPar) %>% 
  summarise(n=length(FemPar)) %>% 
  ungroup()
summary_dads=pheno_data %>% 
  group_by(MalePar) %>% 
  summarise(n=length(MalePar)) %>% 
  ungroup()

family_plot=ggplot(summary_families, aes(x = Family, y = n)) +
  geom_bar(stat = "identity", color = "purple4", fill = "mediumpurple1") +
  labs(title = "Number per of Trees per Family", y = "count", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
family_plot

female_plot=ggplot(summary_moms, aes(x = FemPar, y = n)) +
  geom_bar(stat = "identity", color = "magenta4", fill = "pink") +
  labs(title = "Progeny from Each Female Parent", y = "count", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
female_plot

male_plot=ggplot(summary_dads, aes(x = MalePar, y = n)) +
  geom_bar(stat = "identity", color = "steelblue", fill = "lightskyblue1") +
  labs(title = "Progeny from Each Male Parent", y = "count", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
male_plot

###make a figure
family_plot/(female_plot+male_plot)
```

## Family Values <3
Now, I want to show how the phenotyping values vary by family.

### pHK
```{r}
family_pHK=pheno_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(n=sum(!is.na(pHK)), 
                   mean=mean(pHK), 
                   sd=sd(pHK), 
                   se=sd/sqrt(n), 
                   lower.CI=mean-(1.96*se), 
                   upper.CI=mean+(1.96*se)) %>% 
  ungroup()

family_pHK_plot=ggplot(family_pHK, aes(x = Family, y = mean)) +
  geom_line(size = 1.5, color = "magenta4") +
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI),
                width = 0.2, color = "magenta4") +
  geom_point(size = 3, shape = 21, fill = "pink", color = "magenta4") +
  labs(x = "Family", y = "pHK") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
family_pHK_plot

ggplot(pheno_data, aes(x = pHK)) +
  geom_histogram(color = "magenta4", fill = "pink", bins = 10) +
  facet_wrap(~ Family) +
  labs(y = "Count", x = "pHK") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

### pL34
```{r}
family_pL34=pheno_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(n=sum(!is.na(pL34)), 
                   mean=mean(pL34), 
                   sd=sd(pL34), 
                   se=sd/sqrt(n), 
                   lower.CI=mean-(1.96*se), 
                   upper.CI=mean+(1.96*se)) %>% 
  ungroup()

family_pL43_plot=ggplot(family_pL34, aes(x = Family, y = mean)) +
  geom_line(size = 1.5, color = "purple4") +
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI),
                width = 0.2, color = "purple4") +
  geom_point(size = 3, shape = 21, fill = "mediumpurple1", color = "purple4") +
  labs(x = "Family", y = "pL34") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
family_pL43_plot

ggplot(pheno_data, aes(x = pL34)) +
  geom_histogram(color = "purple4", fill = "mediumpurple1", bins = 10) +
  facet_wrap(~ Family) +
  labs(y = "Count", x = "pL34") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

### RlzDose
```{r}
family_RlzDose=pheno_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(n=sum(!is.na(RlzDose)), 
                   mean=mean(RlzDose), 
                   sd=sd(RlzDose), 
                   se=sd/sqrt(n), 
                   lower.CI=mean-(1.96*se), 
                   upper.CI=mean+(1.96*se)) %>% 
  ungroup()

family_RlzDose_plot=ggplot(family_RlzDose, aes(x = Family, y = mean)) +
  geom_line(size = 1.5, color = "aquamarine3") +
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI),
                width = 0.2, color = "aquamarine3") +
  geom_point(size = 3, shape = 21, fill = "aquamarine", color = "aquamarine3") +
  labs(x = "Family", y = "RlzDose") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
family_RlzDose_plot
```

### LgLarv
```{r}
family_LgLarv=pheno_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(n=sum(!is.na(LgLarv)), 
                   mean=mean(LgLarv), 
                   sd=sd(LgLarv), 
                   se=sd/sqrt(n), 
                   lower.CI=mean-(1.96*se), 
                   upper.CI=mean+(1.96*se)) %>% 
  ungroup()

family_LgLarv_plot=ggplot(family_LgLarv, aes(x = Family, y = mean)) +
  geom_line(size = 1.5, color = "yellowgreen") +
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI),
                width = 0.2, color = "yellowgreen") +
  geom_point(size = 3, shape = 21, fill = "darkolivegreen1", color = "yellowgreen") +
  labs(x = "Family", y = "LgLarv") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
family_LgLarv_plot
```

### figure
```{r}
family_pHK_plot / family_pL43_plot
```

## Species Distributions
Now, I want to take a look at the absolute disaster that is the species groupings.

### Add the species info
I created a column called "Species-ish" which lists if it is a hybrid and what kind. There are a lot of combos here.
```{r}
###import the data
species=read.csv("USFS_Tails_Metadata_2023.xlsx - Parentage_Species_IDs(Marne'sVersion).csv", header = T)

###add of the species info
species_pheno_data=left_join(pheno_data, species, by = "Family")
```

### pHK
```{r}
species_pHK=species_pheno_data %>% 
  group_by(Species.ish) %>% 
  dplyr::summarise(n=sum(!is.na(pHK)), 
                   mean=mean(pHK), 
                   sd=sd(pHK), 
                   se=sd/sqrt(n), 
                   lower.CI=mean-(1.96*se), 
                   upper.CI=mean+(1.96*se)) %>% 
  ungroup()

species_pHK_plot=ggplot(species_pHK, aes(x = Species.ish, y = mean)) +
  geom_line(size = 1.5, color = "magenta4") +
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI),
                width = 0.2, color = "magenta4") +
  geom_point(size = 3, shape = 21, fill = "pink", color = "magenta4") +
  labs(x = "Species", y = "pHK") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
species_pHK_plot

species_pHK_facet=ggplot(species_pheno_data, aes(x = pHK)) +
  geom_histogram(color = "magenta4", fill = "pink", bins = 10) +
  facet_wrap(~ Species.ish) +
  labs(y = "Count", x = "pHK") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
species_pHK_facet
```

### pL34
```{r}
species_pL34=species_pheno_data %>% 
  group_by(Species.ish) %>% 
  dplyr::summarise(n=sum(!is.na(pL34)), 
                   mean=mean(pL34), 
                   sd=sd(pL34), 
                   se=sd/sqrt(n), 
                   lower.CI=mean-(1.96*se), 
                   upper.CI=mean+(1.96*se)) %>% 
  ungroup()

species_pL43_plot=ggplot(species_pL34, aes(x = Species.ish, y = mean)) +
  geom_line(size = 1.5, color = "purple4") +
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI),
                width = 0.2, color = "purple4") +
  geom_point(size = 3, shape = 21, fill = "mediumpurple1", color = "purple4") +
  labs(x = "Species", y = "pL34") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
species_pL43_plot

species_pL34_facet=ggplot(species_pheno_data, aes(x = pL34)) +
  geom_histogram(color = "purple4", fill = "mediumpurple1", bins = 10) +
  facet_wrap(~ Species.ish) +
  labs(y = "Count", x = "pL34") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
species_pL34_facet
```

### making figures
```{r}
(species_pHK_plot + species_pL43_plot) / (species_pHK_facet + species_pL34_facet)
```
That's not a great figure but they look interesting on their own.

## What do we need?
So we have a fair number of individuals that don't have phenotype or species data. I need to identify what individuals were missing phenotype data for, and which families don't have species data.

### No phenotype data
```{r}
###join the dataframes
pheno_data_withNA=left_join(inds_seq, pheno_data, by = "IndShortName")

###remove the inds that we have sequence data for but no phenotype data
missing_pheno=pheno_data_withNA %>% 
  filter(is.na(pHK)) %>% 
  dplyr::select(IndShortName)

###save the csv
write.csv(missing_pheno, "C:/Users/mquig/Desktop/ash/GWAS_tails/phenotype/missing_phenotypes.csv")
```

### no species data
```{r}
missing_species=species_pheno_data %>% 
  filter(is.na(Species.ish)) %>% 
  dplyr::select(IndShortName, Family)
```
</details>

<details>
  <summary>Phenotype Visuals</summary>

<img width="850" height="768" alt="linear_models" src="https://github.com/user-attachments/assets/cf24bbe4-46e5-4df6-ab8e-6995e886a4d9" />

As evidenced by the linear models, pHK and pL34 have the most significant relationships. Then, I checked to see how these values are distributed across families.

<img height="500" alt="mean and 95% CI plots of pHK and pL34 per family" src="https://github.com/user-attachments/assets/0a676705-e17d-4490-af9e-e3c1e284a580" />

Some families clearly score better than others. We want a high pHK and a low pL34. This means the tree is killing off the larvae before they reach instars 3 and 4.
</details>

# Attempt 1: 369 samples (started 7/27/2026)
For my first GWAS attempt, I am going to run it on the 369 individuals that we have both genotype and phenotype data on. This includes greenish and whiteish ash, and both sexes. I am putting all of the associated script in a parent folder called **attempt1_7-27-2026**.

## Step 2: SNP Calling
<details>
   <summary>1) First, I moved the files with the slurm script *00.move_files.sh*.</summary>

   ```bash
   #!/bin/bash
   #SBATCH --job-name=move_analysis
   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=15
   #SBATCH --mem=32G
   #SBATCH -A ACF-UTK0011
   #SBATCH --partition=short
   #SBATCH --qos=short
   #SBATCH --output=logs/cp_fqgz_%j.out
   #SBATCH --error=logs/cp_fqgz_%j.err
   #SBATCH --time=2:00:00
   #SBATCH --mail-type=BEGIN,END,FAIL
   #SBATCH --mail-user=mquigg1@vols.utk.edu

   export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

   # ─── USER CONFIGURATION ──────────────────────────────────────────────────────
   SOURCE_DIR="/lustre/isaac24/proj/UTK0032/TSIP_ash/NRS_data_master"
   DEST_PATH="/lustre/isaac24/scratch/mquigg1/tails_GWAS/00.input_data"
   PARALLEL_JOBS=15
   # ─────────────────────────────────────────────────────────────────────────────

   # Activate conda environment properly in a non-interactive shell
   source "$(conda info --base)/etc/profile.d/conda.sh"
   conda activate basics

   # Create log dir if it doesn't exist
   mkdir -p logs

   echo "Job started: $(date)"
   echo "Searching for .fq.gz files under: ${SOURCE_DIR}"

   # Find all .fq.gz files recursively (following symlinks with -L),
   # then copy each file in parallel, dereferencing symlinks with rsync -L
   find -L "${SOURCE_DIR}" -type f -name "*.fq.gz" | \
       parallel --jobs "${PARALLEL_JOBS}" \
                --plain \
                --joblog logs/parallel_transfer_log_${SLURM_JOB_ID}.txt \
                --resume \
       "rsync -aL {} ${DEST_PATH}/"

   EXIT_CODE=$?

   if [ ${EXIT_CODE} -eq 0 ]; then
       echo "All copies completed successfully."
   else
       echo "One or more copies failed. Check logs/parallel_transfer_log_${SLURM_JOB_ID}.txt for details."
   fi

   echo "Job finished: $(date)"
   exit ${EXIT_CODE}
   ```

   </details>

  <details>
   <summary>2) I used Claude to write a script that takes the sampleID from a csv file and moves it into a folder. In this case, I made a csv file with all of the individuals that we have phenotype data for and directed them into a new folder. This is script *00.move_files_csv.sh*. This does not use slurm, so execute it [bash organize_fastq.sh ./fastq_files samples.csv ./matched_samples].</summary>

   ```bash
  #!/bin/bash

# Usage: bash organize_fastq.sh <fastq_dir> <sample_csv> <output_dir>
# Example: bash organize_fastq.sh ./fastq_files samples.csv ./matched_samples

FASTQ_DIR="${1}"
SAMPLE_CSV="${2}"
OUTPUT_DIR="${3:-./matched_samples}"

# --- Validate inputs ---
if [[ -z "$FASTQ_DIR" || -z "$SAMPLE_CSV" ]]; then
    echo "Usage: bash organize_fastq.sh <fastq_dir> <sample_csv> [output_dir]"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

if [[ ! -f "$SAMPLE_CSV" ]]; then
    echo "ERROR: Sample CSV not found: $SAMPLE_CSV"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "FASTQ directory : $FASTQ_DIR"
echo "Sample CSV      : $SAMPLE_CSV"
echo "Output directory: $OUTPUT_DIR"
echo "----------------------------------------"

matched=0
missing=0

# Read each sample ID from the CSV (handles single-column or multi-column CSVs)
# Skips blank lines and strips carriage returns / surrounding whitespace
while IFS=',' read -r sample_id _rest; do
    # Clean up the sample ID
    sample_id=$(echo "$sample_id" | tr -d '\r' | xargs)

    # Skip empty lines or header-like lines if needed
    [[ -z "$sample_id" ]] && continue

    r1="${FASTQ_DIR}/${sample_id}_R1.fq.gz"
    r2="${FASTQ_DIR}/${sample_id}_R2.fq.gz"

    found=0

    if [[ -f "$r1" ]]; then
        cp "$r1" "$OUTPUT_DIR/"
        echo "  [COPIED] $(basename "$r1")"
        ((found++))
    fi

    if [[ -f "$r2" ]]; then
        cp "$r2" "$OUTPUT_DIR/"
        echo "  [COPIED] $(basename "$r2")"
        ((found++))
    fi

    if [[ $found -gt 0 ]]; then
        ((matched++))
    else
        echo "  [MISSING] No files found for sample: $sample_id"
        ((missing++))
    fi

done < "$SAMPLE_CSV"

echo "----------------------------------------"
echo "Done. Matched: $matched sample(s) | Missing: $missing sample(s)"
echo "Files copied to: $OUTPUT_DIR"
   ```

   </details>

  <details>
   <summary>3) Next, I made a metadata table. Vary_cool is vary_particular about the metadata file so I copied some script from Zane. It is in *01.create.metadata*. That metadata file needs to be in the parent GWAS folder so it can be accessed by nextflow.</summary>

   ```bash
# Initialize sample names
ls -1 yes_phenotype_data/*_R1.fq.gz | sed 's/yes_phenotype_data\///g' | sed 's/_R1.fq.gz//g' > sample_names.txt.temp

# Add header
echo sample,r1,r2 > metadata.csv

# add file names
awk -v DIR="$PWD/yes_phenotype_data/" 'BEGIN {OFS=","} {print $1, DIR $1 "_R1.fq.gz", DIR $1 "_R2.fq.gz"}' sample_names.txt.temp >> metadata.csv

# Clean Up
rm sample_names.txt.temp
   ```

   </details>

<details>
   <summary>4) Finally, I ran the script. It took forever to make it to the top of the queue. The script is *02.run_vary_cool.sh*. Update: I switched to long-bigmem and it started running on 7/28/2026 @ 10:50am.</summary>

   ```bash
#!/bin/bash
#SBATCH -J SNP_call_attempt2
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6-00:00:00
#SBATCH --error=logs/job.SNP_call.e%J
#SBATCH --output=logs/job.SNP_call.o%J
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mquigg1@vols.utk.edu

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate vary_cool

export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOG=false

nextflow /lustre/isaac24/scratch/mquigg1/tails_GWAS/00.vary_cool/vary_cool/main.nf \
    --publish_dir /lustre/isaac24/scratch/mquigg1/tails_GWAS/01.vary_cool_output \
    --input /lustre/isaac24/scratch/mquigg1/tails_GWAS/metadata.csv \
    --genome /lustre/isaac24/scratch/mquigg1/tails_GWAS/00.references/pe57_v-T-B.H.C.C.A.A.FINAL.hap1.fasta \
    --skip_qc false \
    --skip_trim false \
    --skip_mark_dupe false \
    --aligner bwa_mem \
    --ploidy 2 \
    --chunks 12 \
    --caller bcftools \
    --bp_intervals 10000000 \
    -profile slurm,custom \
    -resume
   ```

   </details>


## Run GWAS
Zane recommended two different options for analysis.
  1) [vcf2gwas with GEMMA](https://github.com/frankvogt/vcf2gwas)
  2) [GAPIT](https://github.com/jiabowang/GAPIT)
I'm going to start with option 1) vcf2gwas.

<details>
  <summary>Options</summary>
  -v / --vcf
Specify genotype .vcf or .vcf.gz file (required).

-pf / --pfile
Specify phenotype file.

-p / --pheno
Specify phenotypes used for analysis:
Type the phenotype name
OR
'1' selects first phenotype from phenotype file (second column), '2' the second phenotype (third column) and so on.

-ap / --allphenotypes
All phenotypes in the phenotype file will be used.

-cf / --cfile
Type 'PCA' to extract principal components from the VCF file
OR
Specify covariate file.

-c / --covar
If 'PCA' selected for the -cf / --cfile option, set the amount of PCs used for the analysis
Else:
Specify covariates used for analysis:
Type the covariate name
OR
'1' selects first covariate from covariate file (second column), '2' the second covariate (third column) and so on.

-ac / --allcovariates
All covariates in the covariate file will be used.

-chr/ --chromosome
Specify chromosomes for analysis.
By default, all chromosomes will be analyzed.
Input value has to be in the same format as the CHROM value in the VCF file

-gf / --genefile
Specify gene file.

-gt / --genethresh
Set a gene distance threshold (in bp) when comparing genes to SNPs from GEMMA results.
Only SNPs with distances below threshold will be considered for comparison of each gene.

-k / --relmatrix
Specify relatedness matrix file.

-o/ --output
Change the output directory.
Default is the current working directory.

GEMMA affiliated options:
-lm {1,2,3,4}
Association Tests with a Linear Model.
optional: specify which frequentist test to use (default: 1)
1: performs Wald test
2: performs likelihood ratio test
3: performs score test
4: performs all three tests

-gk {1,2}
Estimate Relatedness Matrix from genotypes.
optional: specify which relatedness matrix to estimate (default: 1)
1: calculates the centered relatedness matrix
2: calculates the standardized relatedness matrix

-eigen
Perform Eigen-Decomposition of the Relatedness Matrix.

-lmm {1,2,3,4}
Association Tests with Univariate Linear Mixed Models.
optional: specify which frequentist test to use (default: 1)
1: performs Wald test
2: performs likelihood ratio test
3: performs score test
4: performs all three tests
To perform Association Tests with Multivariate Linear Mixed Models, set '-multi' option

-bslmm {1,2,3}
Fit a Bayesian Sparse Linear Mixed Model
optional: specify which model to fit (default: 1)
1: fits a standard linear BSLMM
2: fits a ridge regression/GBLUP
3: fits a probit BSLMM

-m / --multi
performs multivariate linear mixed model analysis with specified phenotypes
only active in combination with '-lmm' option

-w / --burn
specify burn-in steps when using BSLMM model.
Default value: 100,000

-s / --sampling
specify sampling steps when using BSLMM model.
Default value: 1,000,000

-smax / --snpmax
specify maximum value for 'gamma' when using BSLMM model.
Default value: 300

Miscellaneous options:
-M / --memory
set memory usage (in MB)
if not specified, half of total memory will be used

-T / --threads
set core usage
if not specified, all available logical cores minus 1 will be used

-q / --minaf
minimum minor allele frequency (MAF) of sites to be used (default: 0.01)
input value needs to be a value between 0.0 and 1.0

-ts / --topsnp
number of top SNPs of each phenotype to be summarized (default: 15)
after analysis the specified amount of top SNPs from each phenotype will be considered

-P / --PCA
perform PCA on phenotypes and use resulting PCs as phenotypes for GEMMA analysis
optional: set amount of PCs to be calculated (default: 2)
recommended amount of PCs: 2 - 10

-U / --UMAP
perform UMAP on phenotypes and use resulting embeddings as phenotypes for GEMMA analysis
optional: set amount of embeddings to be calculated (default: 2)
recommended amount of embeddings: 1 - 5

-um / --umapmetric
choose the metric for UMAP to use to compute the distances in high dimensional space
Default: euclidean
Available metrics: euclidean, manhattan, braycurtis, cosine, hamming, jaccard, hellinger

-t / --transform
transform the input phenotype file
applies the selected metric across rows
Default: wisconsin
Available metrics: total, max, normalize, range, standardize, hellinger, log, logp1, pa, wisconsin

-asc / --ascovariate Use dimensionality reduction of phenotype file via UMAP or PCA as covariates
Only works in conjunction with -U / --UMAP or -P / --PCA

-KC / --kcpca Kinship calculation via principal component analysis instead of GEMMA's internal method
optional: r-squared threshold for LD pruning (default: 0.5)

-sv / --sigval
set value where to draw significant line in manhattan plot
represents -log10(1e-).
Default: Bonferroni corrected with total amount of SNPs used for analysis.
set to '0' to disable line

-nl / --nolabel
remove the SNP labels in the manhattan plot
reduces runtime if analysis results in many significant SNPs

-nq / --noqc
deactivate Quality Control plots
reduces runtime

-np / --noplot
deactivate Manhattan and QQ-plots
reduces runtime

-fs/ --fontsize
set the fontsize of plots.
Default value: 26

-sd / --seed
perform UMAP with random seed
reduces reproducibility

-r / --retain
keep all temporary intermediate files
e.g. subsetted and filtered VCF and .csv files
</details>
