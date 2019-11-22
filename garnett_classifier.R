suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(garnett))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(argparse))

#################################
Seurat2Monocle<-function(seurat,slot="counts"){
  #Load Seurat object
  require(monocle)
  if(class(seurat)[[1]]=="Seurat"){
    seurat_object<-seurat
  }else{
    seurat_object <- readRDS(seurat)
  }

  #Extract data, phenotype data, and feature data from the SeuratObject
  #data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')
  #data<-as.sparse(seurat_object@assays$RNA@data)
  data<-GetAssayData(seurat_object,slot)
  pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)

  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

  #Construct monocle cds
  print("Create monocle object")
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.1,
                                expressionFamily = negbinomial.size())
  #print("save monocle object")
  #saveRDS(monocle_cds,file=monocle.path)
  return(monocle_cds)
}
################################

print("--------Configure Parameters--------")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--dataset",
                    type="character",
                    default="",
                    help="the dataset  to be used")


parser$add_argument("--classifier",
                    type="character",
                    default="",
                    help="the dataset  to be used")


parser$add_argument("--markers_file",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--train",
                    dest="train",
                    action="store_true")

args<-parser$parse_args()
#  Load object
model.dir<-file.path(args$dataset,"model","seurat.rds")
plot.dir<-file.path(args$dataset,"plot")

##########################
print("###  Load Monocle Object")
seurat<-readRDS(model.dir)
print("### Convert ...")
cds<-Seurat2Monocle(seurat,"counts")

print("### estimateSizeFactors")
cds<-estimateSizeFactors(cds)
classifier<-readRDS(args$classifier)
markers_file<-args$markers_file

###########################

if(args$train){
	print("*** Training ***")
	classifier <- train_cell_classifier(cds = cds,
                                         marker_file = markers_file,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
}

print("*** Use pretrained classifier to predict cell types ***")
cds <- classify_cells(cds,classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
			   verbose=TRUE,
                           cds_gene_id_type = "SYMBOL")

print("### Saving")
saveRDS(cds,file.path(args$dataset,"model","garnnet_cds.rds"))
seurat@meta.data$garnett_cluster<-pData(cds)$cell_type
saveRDS(seurat,model.dir)
print("### Done")
