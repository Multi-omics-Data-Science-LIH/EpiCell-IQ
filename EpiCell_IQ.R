
source("housman_test_cfDNA.R")

#Read in

library(methrix)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(RColorBrewer)

##download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185307

hg38 <- methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

# Step 1: Different fragment size

# Fragment size 160-180bp
meth_lung_160_180 <- methrix::read_bedgraphs(
  files = bed_files_160_180,
  ref_cpgs = hg38,
  stranded = TRUE,
  collapse_strands = TRUE,
  chr_idx = 1, start_idx = 2,
  end_idx = 3, U_idx = 12, M_idx = 13
)

# Fragment size 240-400bp
meth_lung_240_400 <- methrix::read_bedgraphs(
  files = bed_files_240_400,
  ref_cpgs = hg38,
  stranded = TRUE,
  collapse_strands = TRUE,
  chr_idx = 1, start_idx = 2,
  end_idx = 3, U_idx = 12, M_idx = 13
)

# Fragment size 400bp
meth_lung_400 <- methrix::read_bedgraphs(
  files = bed_files_400,
  ref_cpgs = hg38,
  stranded = TRUE,
  collapse_strands = TRUE,
  chr_idx = 1, start_idx = 2,
  end_idx = 3, U_idx = 12, M_idx = 13
)

## remove_uncovered
filt_160-180 <- remove_uncovered(meth_lung_160_180)
filt_240-400 <- remove_uncovered(meth_lung_240_400)
filt_400 <- remove_uncovered(meth_lung_400)

#methylation atlas from:
#Moss, J., Magenheim, J., Neiman, D. et al.
#Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease.
#Nat Commun 9, 5068 (2018). https://doi.org/10.1038/s41467-018-07466-6


#https://github.com/nloyfer/meth_atlas/blob/master/full_atlas.csv.gz
atlas <- read.csv("reference_meth/atlas.csv", row.names=1)
atlas <- atlas[-1,]

### 450k manifest
### downloaded from: https://zwdzwd.github.io/InfiniumAnnotation
HM450.hg38.manifest <- read.delim("reference_meth/HM450.hg38.manifest.tsv.gz")
rownames(HM450.hg38.manifest) <- HM450.hg38.manifest$Probe_ID
HM450.hg38.manifest$CpG_beg <- HM450.hg38.manifest$CpG_beg + 1
colnames(HM450.hg38.manifest)[1:3]  <- c("chr", "start", "end")
HM450.hg38.manifest <- HM450.hg38.manifest[HM450.hg38.manifest$chr %in% paste0("chr", 1:22),]
common_cpg <- intersect(rownames(HM450.hg38.manifest), rownames(atlas))

HM450.hg38.manifest <- HM450.hg38.manifest[common_cpg,]
HM450 <- makeGRangesFromDataFrame(HM450.hg38.manifest)
HM450 <- HM450[!duplicated(HM450),]

common_cpg <- intersect(rownames(HM450.hg38.manifest), rownames(atlas))
atlas <- atlas[common_cpg,]

atlas_450k <- methrix:::create_methrix(beta_mat = atlas, cov_mat = ceiling(atlas), cpg_loci = HM450.hg38.manifest[common_cpg, c("chr", "start", "end")],
                                       is_hdf5 = FALSE, genome_name = "hg38", col_data = data.frame(samples=colnames(atlas), row.names=colnames(atlas)), h5_dir = NULL,
                                       ref_cpg_dt = hg38$cpgs, chrom_sizes = hg38$contig_lens, desc = NULL)
atlas_450k <- atlas_450k[names(HM450),]

## Pivot epithelial

pivot_base <- read.csv("reference_meth/epithelial.csv")
pivot_base$destination[pivot_base$source == "Lung cells"] <- "Lung cells"     # Extract Lung cells sample from Epithlial 
pivot_base$destination[pivot_base$destination==""] <- "Other"
pivot <- data.frame(matrix(NA, ncol=length(unique(pivot_base$destination)), nrow = nrow(pivot_base)))
colnames(pivot) <- unique(pivot_base$destination)
rownames(pivot) <- pivot_base$source

for (i in seq_along(1:nrow(pivot))){
  pivot[i,]  <- ifelse(pivot_base$destination[i]==colnames(pivot), 1, 0)}


rownames(pivot) <- gsub(" ", "_", rownames(pivot))
rownames(pivot) <- gsub("-", ".", rownames(pivot))
colnames(pivot) <- gsub(" ", "_", colnames(pivot))


## Ordering the most variable sites
atlas2_450k_most_variable <- list(
  "50k" = methrix::order_by_sd(atlas_450k)[1:50000,],
  "20k" = methrix::order_by_sd(atlas_450k)[1:20000,],
  "5k" = methrix::order_by_sd(atlas_450k)[1:5000,],
  "1k" = methrix::order_by_sd(atlas_450k)[1:1000,]
)

# Initialize results list
results_atlas2 <- list()

# Create a list of parameters for each analysis
analysis_params_atlas2 <- list(
  list(i = "atlas2_180", frag = 180, meth_data = filt_160-180),
  list(i = "atlas2_240_400", frag = "240_400", meth_data = filt_240-400 ),
  list(i = "atlas2_400", frag = 400, meth_data = filt_400)
)

# Loop through each set of parameters for both 100k and 10k most variable CpGs
overlapping_list_atlas2 <- list()
transformed_datat_atlas2 <- list()  # To store the transformed data

for (params in analysis_params_atlas2) {
  # Extract parameters
  i <- params$i
  frag <- params$frag
  meth_data <- params$meth_data
  
  for (var_set in c("50k", "20k" ,"5k", "1k")) {
    # Perform the analysis for each variable set ("50k", "20k" ,"5k" and "1k")
    res_atlas2 <- houseman_cell_type(m = meth_data, reference_meth = atlas2_450k_most_variable[[var_set]], 
                                     included_regions = NULL, pivot = pivot, 
                                     included_cell_types = NULL)
    
    # Store the results for the respective most variable set
    results_atlas2[[i]][[paste0("most_", var_set, "_variable_", frag)]] <- round(res_atlas2[[1]], 3)
    
    # Store the overlapping results
    overlapping_list_atlas2[[i]][[var_set]] <- res_atlas2[[2]]
    
    # Dynamically transform the data and store it
    transformed_datat_atlas2[[paste0("data_", var_set, "_", frag)]] <- t(results_atlas2[[i]][[paste0("most_", var_set, "_variable_", frag)]])
  }
}

# Now the results are stored for most variable CpGs

##Plot heatmaps

library(circlize)
library(ComplexHeatmap)

# Define the color function
col_fun <- colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred"))

# Define the sample type annotations
sample_types <- c(rep("Cancer", 6), rep("Healthy", 7))

# Define the annotations for the heatmaps
col_anno <- HeatmapAnnotation(
  sample = sample_types,
  col = list(sample = c("Cancer" = "red", "Healthy" = "green"))
)

# Extract fragments, ensuring both numeric and underscored fragments like 240_400 are included
fragments <- unique(sub(".*_(\\d+_\\d+|\\d+)$", "\\1", names(transformed_data_atlas_450k)))

# Add 240_400 explicitly if not included
if (!"240_400" %in% fragments) {
  fragments <- c(fragments, "240_400")
}

# Print fragments to verify inclusion of 240_400
print(fragments)

# Define the CpG types
var_sets <- c("50k", "5k")

# Loop through both fragments and CpG types and draw each heatmap separately
for (frag in fragments) {
  for (var_set in var_sets) {
    # Dynamically retrieve the transformed data
    data_key <- paste0("data_", var_set, "_", frag)
    
    # Check if the dataset exists
    if (!is.null(transformed_data_atlas_450k[[data_key]])) {
      cat("Processing:", data_key, "\n")  # Debugging output
      
      # Create the heatmap
      heatmap <- Heatmap(
        transformed_data_atlas_450k[[data_key]],
        cluster_rows = FALSE, cluster_columns = FALSE, 
        top_annotation = col_anno, show_column_names = FALSE,
        row_title = "All sites", 
        column_title = paste("Highly variable sites Atlas 450k (", var_set, " CpGs), fragment size ", frag),
        col = col_fun
      )
      
      # Draw the heatmap separately
      draw(heatmap, heatmap_legend_side = "right")
    } else {
      cat("Data for key", data_key, "not found!\n")  # Debugging output
    }
  }
}









