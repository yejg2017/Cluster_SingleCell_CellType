suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))

#####################
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

######################
print("--------Configure parameters--------")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--dataset",
                    type="character",
                    default="",
                    help="which dataset to use")

args <- parser$parse_args()
###############################
outdir<-dirname(args$dataset)
print(paste0("Loading dataset from : ",args$dataset))
print(paste0("Result will be saved into ",outdir))


print("Step 1 : Loading data")
object<-readRDS(args$dataset)

print("Step 2 : Convert into newcellDataset object")
cds<-Seurat2Monocle(object,"data")

###############################
print("Step 3 : Create a CellTypeHierachy")
cth <- newCellTypeHierarchy()

#cth<-addCellType(cth,"CD45+ Cell",
#		 classify_func=function(x){x["PTPRC",]>0})

cth<-addCellType(cth,"T Cell",
		 classify_func=function(x){x["CD3D",]>0})

cth<-addCellType(cth, "DCs",
                 classify_func=function(x) { x["IL3RA",] > 0 | x["ITGAX",]>0})

cth<-addCellType(cth, "pDCs",
                 classify_func=function(x) { x["GZMB",]>0},
		 parent_cell_type_name="DCs")

cth<-addCellType(cth, "cDC1",
                 classify_func=function(x) { x["CLEC9A",]>0},
                 parent_cell_type_name="DCs")

cth<-addCellType(cth, "pre-DC",
                 classify_func=function(x) { x["AXL",]>0 & x["SIGLEC6",]>0},
                 parent_cell_type_name="DCs")

cth<-addCellType(cth, "CLEC10A+ CLEC4A low DC2",
                 classify_func=function(x) {x["CD1C",]>0 & x["CLEC10A",]>0 & x["CLEC4A",]<quantile(x["CLEC4A",])[2]},
                 parent_cell_type_name="DCs")

cth<-addCellType(cth, "CLEC10A+ CLEC4A Hi DC2",
                 classify_func=function(x) { x["CD1C",]<quantile(x["CD1C",])[2] & x["CLEC10A",]==0 & x["CLEC4A",]>quantile(x["CLEC4A",])[3]},
                 parent_cell_type_name="DCs")

cth <- addCellType(cth, "CD4+ T Cell", 
                   classify_func=function(x) {x["CD4",] > 0}, 
                   parent_cell_type_name = "T Cell")

cth <- addCellType(cth, "CD8+ T Cell", 
                   classify_func=function(x) {
                     x["CD8A",] > 0 | x["CD8B",] > 0
                   }, 
                   parent_cell_type_name = "T Cell")

cth <- addCellType(cth, "B Cell",
                   classify_func=function(x) {x["CD19",] > 0 | x["MS4A1",]>0})

cth <- addCellType(cth, "Plasma Cell",
                   classify_func=function(x) {x["CD38",] > 0},
                   parent_cell_type_name="B Cell")


cth <- addCellType(cth, "Neutrophil",
                   classify_func=function(x) {x["FUT4",]>0})

cth <- addCellType(cth, "NK Cell",
                   classify_func=function(x) {x["NCAM1",]>0})

cth <- addCellType(cth, "NK T Cell",
                   classify_func=function(x) {x["ZNF683",]>0},
                   parent_cell_type_name="T Cell")

cth <- addCellType(cth, "Monocyte",
                   classify_func=function(x) {x["CD14",]>0 | x["FCGR3A",]>0})

cth <- addCellType(cth, "Macrophage",
                   classify_func=function(x) {x["CD163",]>0 | x["LYZ",]>0})



print("Step 4 : Train Classifier")
cds<-classifyCells(cds,cth,frequency_thresh=0.01)

print("Step 5 : Save the update object")
saveRDS(cds,file.path(outdir,"cds.rds"))
print("Done!!!")





