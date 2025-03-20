# Birth-Death Model of Gene Family Evolution and Detecting Whole-Genome Multiplications (WGMs) using Reciprocally Retained Gene Families

## About

Given a species tree (with annotated Whole-Genome Multiplication (WGM) nodes) and a gene family size profile across species, this model estimates the likelihood of the observed gene family sizes based on a duplication/loss rate (λ) and an ancestral gene family size (r) at the root of the tree. The likelihood is maximized as a function of λ, for fixed values of r ranging from 1 to 20, using a cutting-plane optimization method.

Additionally, this repository provides tools for reverse-engineering WGMs using a combination of stochastic birth-death (BD) modeling and a complementary log-likelihood ratio test (cLRT). The cLRT pipeline leverages reciprocally retained (RR) gene families, inferred through BD modeling, as WGM markers to detect or reject WGM events along branches of a phylogenetic tree.

For more information regarding the proposed mathematical models please see my PhD thesis (http://hdl.handle.net/1854/LU-8637576) chapter 2 and chapter 4, 
Chapter2: Reciprocally Retained Genes in the Angiosperm Lineage Show the Hall-marks of Dosage Balance Sensitivity (https://pubmed.ncbi.nlm.nih.gov/29061868/)
Chapter 4: Reverse Engineering of WGMs

The project includes implementations of:
- **BirthDeathModel**: Computes the likelihood of gene family size distributions given duplication/loss rates (BD model).
- **WGM Detection Framework**: Uses marker GFs (output of BD modeling) to detect WGM events based on gene family retention patterns (cLRT pipeline).

This repository is particularly useful for researchers studying plant genome evolution, polyploidy, WGMs, gene retention dynamics and dosage balance effects. 

## Authors
This code was developed by:
- **Setareh Tasdighian** (Ghent University, 9000 Ghent, Belgium)
- **[@stasdigh](https://github.com/stasdigh)** - original developer

For inquiries, contact: **setareh.tasdighian@gmail.com**

## Installation
To install and set up the project, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/stasdigh/ReverseEngWGDs.git
   cd ReverseEngWGDs
   ```
2. Ensure you have **Maven** installed for Java-based components.
3. Build the executable jar file:
   ```bash
   mvn clean install
   ```
   If the build is successful, the jar file with embedded dependencies will be located in the `target` folder.

## Usage
### BirthDeathModel
This tool computes the likelihood of gene family size distributions given a species tree and gene counts.
#### Input Parameters:
1. **Species tree** in Newick format.
2. **WGM event file** listing whole-genome multiplications in the tree.
3. **Gene family counts file**, specifying gene counts for different species.
4. **Gene family index** (integer) indicating the starting index in the counts file.
5. **Root size** (integer) specifying gene count at the root.

#### Example Usage:
```bash
java -jar target/BirthDeathModel.jar species_tree.txt wgm_events.txt gene_counts.txt 0 10
```
#### Output:
The program returns:
```
Gene_Family_ID \t root_size \t optimal_lambda \t log_likelihood \t p-value
```

### WGM Detection
The framework applies BD modeling to gene family size distributions to infer WGM events.
#### Input:
- Species tree with WGM event placements.
- Gene family size distributions across species.
- Reciprocally retained gene families used as WGM markers.

#### Output:
- **Detected WGM events** with likelihood scores.
- **Number of detected WGMs per marker GF**

## Citation
If you use this repository in your research, please cite:
Tasdighian et al., "Using reciprocally retained gene families to detect whole-genome multiplications in plants." Molecular Biology and Evolution (2025).

## License
This project is released under the GNU General Public License (GPL-3.0).
