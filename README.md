# Cluster_SingleCell_CellType
Use the method to auto classify cell type in single cell

## 这个repos主要是讲述了如何在单细胞（10X）自动化的去识别细胞类型
  - 通过提前定义好marker 条件 或者marker文件

### Requirement
 - R( version >3.5),最好通过 [conda](https://anaconda.org/)来配置,*conda install -c r r-base*
 - Seurat,一个目前比较流行的单细胞分析程序包，参考官网[satijalab](https://satijalab.org/seurat)
 - monocle(version 2)，单细胞分析程序包，主要用来做newCellTypeHierarchy方法的
 - monocle3,这个版本的monocle3看基本上跟monocle（version 2）很大区别，具体看[官网](https://cole-trapnell-lab.github.io/monocle3)
 - garnnet,跟monocle和monocle3配套使用，参考[官网]（https://cole-trapnell-lab.github.io/garnett）
 - ggplot2 ...等等这些不一一赘述
