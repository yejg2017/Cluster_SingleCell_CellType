# Cluster_SingleCell_CellType
Use the method to auto classify cell type in single cell

## 这个repos主要是讲述了如何在单细胞（10X）自动化的去识别细胞类型
  - 通过提前定义好marker 条件 或者marker文件

### Requirement
 - R( version >3.5),最好通过 [conda](https://anaconda.org/)来配置,*conda install -c r r-base*
 - Seurat,一个目前比较流行的单细胞分析程序包，参考官网[satijalab](https://satijalab.org/seurat)
 - cellassign，需要提前装tensorflow,参考官网(https://github.com/Irrationone/cellassign)
 - monocle(version 2)，单细胞分析程序包，主要用来做newCellTypeHierarchy方法的
 - monocle3,这个版本的monocle3看基本上跟monocle（version 2）很大区别，具体看[官网](https://cole-trapnell-lab.github.io/monocle3)
 - garnnet,跟monocle和monocle3配套使用，参考[官网]（https://cole-trapnell-lab.github.io/garnett）
 - ggplot2 ...等等这些不一一赘述
 
###  注意
* Seurat 处理好的seurat 对象目前不能直接通过里面的函数来互相转化了，需要自己码函数
``` R
Seurat2Monocle<-function(seurat,slot="counts"){
  # slot: counts or data or scaled.data
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
```

### Garnett
Garnett方法可以兼容monocle2 和monocle3包，但是，monocle2 需要先使用Seurat2Monocle函数转化为monocle对象才行，
而monocle3不需要，直接使用monocle生成的对象来使用garnett里面的方法（garnett针对monocle2 ，
monocle3的安装方法也稍微不一样，具体看官网安装方法）
   - [garnett_classifier](./garnett_classifier)文件夹保存着提前训练好的garnett 分类器以及实例的marker.txt 文件
