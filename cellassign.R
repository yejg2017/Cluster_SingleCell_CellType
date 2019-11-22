suppressPackageStartupMessages(library(cellassign))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

print("--------Configure parameters--------")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--dataset",
                    type="character",
                    default="",
                    help="dataset to be used")

parser$add_argument("--train",
                    dest="train",
                    action="store_true")

args <- parser$parse_args()
############################
model.dir<-file.path(args$dataset,"model","seurat.rds")
outdir<-file.path(args$dataset,"cellassign")

print(paste0("Loading dataset from : ",model.dir))
print(paste0("Save output in : ",outdir))

if(!dir.exists(outdir)){
	dir.create(outdir)
}

##########################
print("Step 1 : Loading Dataset and Convert into SingleCellExperiment")
seurat<-readRDS(model.dir)


print("Step 3 : Create  binary marker by cell type matrix")
marker_list<-list()
marker_list$`rod cell`<-c("RHO","SAG","RBP1","CNGA1","PDE6A")
marker_list$`cone cell`<-c("GNB3","CLUL1","SLC24A2","DST","TTR")
marker_list$`biopolar cell`<-c("TRPM1","GRM6")
marker_list$`astrocytes`<-c("GFAP")
marker_list$`amacrine cell`<-c("GAD1","CALB1","PROX1","TFAP2A","SLC6A9","BARHL2")
marker_list$`muller cell`<-c("SRGN","RDH5","RLBP1")
marker_list$`horizontal cell`<-c("ONECUT1","ONECUT2")
marker_list$`microglia`<-c("TYROBP","HLA-DRA","HLA-DPB1","HLA-DPA1")
marker_list$`ganglion cell`<-c("NEFL","POU4F1","GAP43","SNCG")


filename<-file.path(outdir,"cellassgin.rds")
if(args$train){
	print("Step 2 : Estimate SizeFactors")
	object<-as.SingleCellExperiment(seurat)
        object<-computeSumFactors(object)
        s<-sizeFactors(object)

	marker_mat <- marker_list_to_mat(marker_list)
        marker_use<-intersect(rownames(marker_mat),rownames(object))

	print("Step 4 : Fit the Cellassign Model")
	fit <- cellassign(exprs_obj = object[marker_use,], 
                  marker_gene_info =marker_mat[marker_use,], 
                  s = s, 
		  sce_assay="counts",
                  learning_rate = 0.1, 
                  shrinkage = TRUE,
                  verbose = TRUE)
	print(fit)
        print("Step 5 : Save the fitted model")
        saveRDS(fit,filename)
}else{
	print("Loading the pretrained cellassign model")
	fit<-readRDS(filename)
}
prediction.cells<-celltypes(fit)
seurat@meta.data$prediction.cells<-prediction.cells
saveRDS(seurat,model.dir)
print("Finish!!!")
