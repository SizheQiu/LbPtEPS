# LbPtEPS: EPS biosynthesis in Lactiplantibacillus plantarum HMX2
## Multi-omics analysis: genomics, metabolomics, and proteomics
### Workflow
1. EPS biosynthetic gene cluster: EPS_BGC.ipynb
2. Metabolomic analysis: Metabolomics.ipynb
3. Differential protein expression analysis: DEG_proteomics.ipynb
4. Analysis of proteome resource allocation under acid stress: ResourceAllocation.ipynb
### Data
1. [Whole genome sequence of L. plantarum HMX2](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Genome_HMX2/genome_hmx2.fa)
2. [Protein sequences of L. plantarum HMX2](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Genome_HMX2/protein_hmx2.fa)
3. [Functional annotations of coding genes](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Genome_HMX2/hmx2_annotation.csv)
4. [EPS biosynthetic gene cluster of L. plantarum HMX2](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Genome_HMX2/BGC_EPS.gbk)
5. [Quantitaive proteomic data](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Proteomics/Proteomics_B.xlsx)
6. [Extra-cellular metabolomic data](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Exp_data/Metabolomics_mM.csv)
7. [Intra-cellular metabolomic data](https://github.com/SizheQiu/LbPtEPS/blob/main/data/Exp_data/IntraMetabolomics.csv)
### Dependencies
1. antiSMASH: https://antismash.secondarymetabolites.org/
2. Biopython: https://biopython.org/
3. DNA Features Viewer: https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
4. PyDESeq2: https://github.com/owkin/PyDESeq2.git
5. Seaborn statistical data visualization:https://seaborn.pydata.org/index.html
6. Scikit-learn: https://scikit-learn.org/

## RPCFBA: regulatory proteome constrained flux balance analysis
### Workflow
1. Genome scale metabolic model modification: GSMM_modification.ipynb
2. Parameter estimation (e.g., the enzyme activity of Glycosyltransferases): Parameter_estimation.ipynb
3. Simulation and in-silico perturbation on carbon sources: Simulation.ipynb
### Dependencies
1. COBRApy: https://github.com/opencobra/cobrapy
2. DLTKcat: https://github.com/SizheQiu/DLTKcat.git
3. Scikit-learn: https://scikit-learn.org/

## Citation
Sizhe Qiu, Aidong Yang, Xinyu Yang, Wenlu Li, Hong Zeng, Yanbo Wang. Proteome trade-off between primary and secondary metabolism shapes acid stress induced bacterial exopolysaccharide production. https://doi.org/10.1101/2024.04.19.590233
