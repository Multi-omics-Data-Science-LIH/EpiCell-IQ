# EpiCell-IQ
EpiCell-IQ is a lightweight R-based tool designed for predicting the originating cell composition (OCC) from methylome data. This tool is particularly optimized for sparse methylation data from Nanopore sequencing but is versatile enough to work with array- and sequencing-based references, making it ideal for identifying the tissue of origin in liquid biopsy applications, including cancer detection.
### Features
- **Flexible reference datasets:** It supports both array- and sequencing-based methylation references, ensuring adaptability to different experimental setups.
- **Efficient deconvolution:** Leverages the computational power of the methrix R package to perform deconvolution using the Houseman algorithm, ensuring fast and accurate results even with sparse data.
- **Cell-type composition estimation:** Provides precise estimates of OCC using highly accurate individual models fitted for each sample, based on non-missing CpG sites.
- **Customizable fragment size analysis:** Explores optimal fragment size intervals for cancer detection, leveraging Nanopore sequencing data.

### Usage
The main function of EpiCell-IQ is in the housman_test_cfDNA.R file. This function is based on the method described in the publication:

Houseman, E.A., Accomando, W.P., Koestler, D.C. et al. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics 13, 86 (2012). https://doi.org/10.1186/1471-2105-13-86

An example use of the tool with cfDNA data is provided in the **EpiCell_IQ.R** script. You can modify this script to suit your dataset.

### Citation
Houseman, E.A., Accomando, W.P., Koestler, D.C. et al. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics 13, 86 (2012). https://doi.org/10.1186/1471-2105-13-86
