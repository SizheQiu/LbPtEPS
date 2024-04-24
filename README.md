# LbPtEPS: EPS biosynthesis in Lactiplantibacillus plantarum HMX2
## Multi-omics analysis: genomics, metabolomics, and proteomics
### Workflow
1. EPS biosynthetic gene cluster: EPS_BGC.ipynb
2. Metabolomic analysis: Metabolomics.ipynb
3. Differential protein expression analysis: DEG_proteomics.ipynb
4. Analysis of proteome resource allocation under acid stress: ResourceAllocation.ipynb
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
