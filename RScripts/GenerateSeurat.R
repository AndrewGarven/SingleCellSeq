options(repos = c(CRAN = "https://cran.r-project.org"))
install.packages('Seurat', lib = '/home/garvena/projects/def-dmberman/garvena/R/library')
.libPaths("~/R/library")
library(Seurat)

# List of sample names and their corresponding data directories
sample_list <- list(
  "SRR12539462" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR12539462/outs/filtered_feature_bc_matrix",
  "SRR12539463" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR12539463/outs/filtered_feature_bc_matrix",
  "SRR12603780" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR12603780/outs/filtered_feature_bc_matrix",
  "SRR14615558" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR14615558/outs/filtered_feature_bc_matrix",
  "SRR9897621" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR9897621/outs/filtered_feature_bc_matrix",
  "SRR9897623" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR9897623/outs/filtered_feature_bc_matrix",
  "SRR9897624" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR9897624/outs/filtered_feature_bc_matrix",
  "SRR9897625" = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/cellrangerlai/SRR9897625/outs/filtered_feature_bc_matrix"
)

# Function to read and create Seurat objects for each sample
read_and_create_seurat <- function(sample_name, data_dir) {
  # Read the data using Read10X
  data <- Read10X(data.dir = data_dir)
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, min.features = 100, project = sample_name)
  
  return(seurat_obj)
}

# Create a list of Seurat objects for each sample
seurat_list <- lapply(names(sample_list), function(sample_name) {
  read_and_create_seurat(sample_name, sample_list[[sample_name]])
})

# Merge all Seurat objects into one
merged_seurat <- merge(x = seurat_list, add.cell.id = names(sample_list))

saveRDS(seurat_obj, file = "/home/garvena/projects/def-dmberman/garvena/singlecellseq/data/LaiSeurat.rds")
