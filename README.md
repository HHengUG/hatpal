# hatpal 🧢

**hatpal** is an R package for alternative polyadenylation (APA) identification & analysis using 3’ scRNA-seq (10x etc.). 

By clustering the possible polyA sites, hatpal can find high-confidence APA clusters without GTF file, and create a count matrix for downstream analysis after being annotated.



## Before the start 📰

### NEW VERSION (hatpal2: https://github.com/HHengUG/hatpal2)
### hatpal2 provides more precise and rich results！

### Dependencies

Softwares: [**samtools**](http://www.htslib.org/download/) [>= 1.9] and [**R**](https://www.r-project.org/) [>= 3.6.1]

R packages: [**data.table**](https://cran.r-project.org/web/packages/data.table/index.html) [1.13.2], [**GenomicAlignments**](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html) [1.24.0], [**kpeaks**](https://cran.r-project.org/web/packages/kpeaks/index.html) [1.1.0] and **[Ckmeans.1d.dp](https://cran.r-project.org/web/packages/Ckmeans.1d.dp/index.html)** [4.3.3]

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
if (!require("data.table")) install.packages("data.table")
if (!require("kpeaks")) install.packages("kpeaks")
if (!require("Ckmeans.1d.dp")) install.packages("Ckmeans.1d.dp")
if (!require("GenomicAlignments")) BiocManager::install("GenomicAlignments")
```



### Installation

The latest hatpal release can be downloaded [here](https://github.com/HHengUG/hatpal/releases).

hatpal can be installed as followed in linux:

```shell
R CMD INSTALL hatpal_0.9.1.tar.gz
```

and in windows:

```R
install.packages("hatpal_0.9.1.zip", repos = NULL, type = "win.binary")
```

hatpal can also be installed from github through [devtools](https://github.com/r-lib/devtools):

```R
library(devtools)
devtools::install_github('HHengUG/hatpal')
```



### Input preparation

- An **indexed BAM** file generated from 3’ scRNA-seq for preprocessing

- A **GTF** file if you need annotation and count matrix for downstream analysis in Seurat

 To find the mapped soft-clipped reads,  run the preprocessing of BAM file on Linux via **samtools**.

```shell
samtools view -F 16 input/input.sorted.bam | awk '{if ($6~/([0-9][0-9]|[3-9])S$/ ) print }' > input/p_3S_test.sam &

samtools view -f 16 input/input.sorted.bam | awk '{if ($6~/^([0-9][0-9]|[3-9])S/ ) print }' > input/n_3S_test.sam &

samtools view -H input/input.sorted.bam > input/SAM_header &

cat input/SAM_header input/p_3S_test.sam input/n_3S_test.sam | samtools view -b > input/3S_input.bam &

samtools sort input/3S_test.bam > input/3S_input.sorted.bam &

samtools index -b input/3S_input.sorted.bam &
```

The output is an indexed BAM file (input/3S_input.sorted.bam) containing soft-clipped reads needed in hatpal.



## Run hatpal 🚴

hatpal includes five steps and functions.

You can simply run the script provided [here](https://github.com/HHengUG/hatpal_example), or use the following functions:

Before you run these functions, you should:

```R
library(hatpal)
```

to library this package

### ExtractChrinfo

*ExtractChrinfo()* extracts the chromosomes information from the indexed BAM file, including name and length.

The input should be a BAM file with a .bai index in the same directory.

```R
chr_info <- ExtractChrinfo("input/input.sorted.bam")
```

The output is a data.table object in R.

### WriteEndlist

*WriteEndlist()* generates two-dimensional coordinate files of each strand for clustering.

Two BAM files and an output directory are needed, and .bai index in the same directory is also required. The next parameter is the output directory, and the last parameter is a data.table containing chromosomes information extracted by *ExtractChrinfo()*.

Use the original indexed BAM file as the first parameter, and the soft-clipped indexed BAM as the second parameter.

```R
WriteEndlist("input/input.sorted.bam", "input/3S_test.sorted.bam", "output/", 
             ExtractChrinfo("input/input.sorted.bam"))
```

The outputs are four files in two directories (negative_strand/ and positive_strand/), and an extra file containing the list of cell barcodes named "output/_CBlist.txt".

### ClusterAPA

*ClusterAPA()* clusters the data in two two-dimensional coordinate files (RC file and SA file) to find the APA clusters.

Every strand has two two-dimensional coordinate files (rc.prebed and sa.prebed) in two directories (negative_strand/ and positive_strand/). These four input files  are  generated by *WriteEndlist()*. The first parameter is the RC file, and the second is the SA file. The last parameter is the output directory. So you should run this function twice for each strand:

```R
ClusterAPA("output/negative_strand/rc.prebed", "output/negative_strand/sa.prebed", 
           "output/negative_strand/")
ClusterAPA("output/positive_strand/rc.prebed", "output/positive_strand/sa.prebed", 
           "output/positive_strand/")
```

Or if you only want high-confidence APA clusters, just replace RC file with SA file:

```R
ClusterAPA("output/negative_strand/sa.prebed", "output/negative_strand/sa.prebed", 
           "output/negative_strand/")
ClusterAPA("output/positive_strand/sa.prebed", "output/positive_strand/sa.prebed", 
           "output/positive_strand/")
```

The outputs are two files containing APA clusters (APA_clusters.out).

### AnnoAPAc

*AnnoAPAc()* annotates the APA clusters from a GTF file.

The first two inputs are  one APA cluster file and one GTF file. The third parameter determines the strand. If chromosome names of the GTF file do not have "chr" in them, the fourth parameter (add_chr = TRUE) will add string "chr" before them. You should also run *AnnoAPAc()* twice for each strand.

```R
AnnoAPAc("output/positive_strand/APA_clusters.out", "ref/hg_19.gtf", "+")
AnnoAPAc("output/negative_strand/APA_clusters.out", "ref/hg_19.gtf", "-")
```

When the chromosome names of the GTF do not have "chr" in it (i.e., 1 2 3), but that of the APA_clusters.out or chr_info have (i.e., chr1 chr2 chr3), run:

```R
AnnoAPAc("output/positive_strand/APA_clusters.out", "ref/hg_19.gtf", "+", add_chr = TRUE)
AnnoAPAc("output/negative_strand/APA_clusters.out", "ref/hg_19.gtf", "-", add_chr = TRUE)
```

The output is a file where each annotation of clusters is described by eight columns. The first common three columns are chr, start and end, and the fourth and fifth columns are annotated gene and type, if there is no annotation for the clusters, then the last three columns will be the gene, type and distance of the closest annotation on the 5' end.

### CountCB

*CountCB()* creates a APA_cluster-Cell_Barcode count matrix which can be used in Seurat analysis.

The first two parameters are annotated APA cluster files of positive strand and negative strand. The third parameter is a BAM file with a .bai index in the same directory. The fourth parameter is a CB list, you can use white list of cell barcode from cellranger, or  use the "output/_CBlist.txt" generated by *WriteEndlist()*. The fifth parameter is a data.table containing chromosomes information extracted by *ExtractChrinfo(),* and the last parameter is the output directory.

```R
CountCB("output/positive_strand/APA_clusters.out.anno", 
        "output/negative_strand/APA_clusters.out.anno", 
        "input/test.sorted.bam","output/_CBlist.txt", 
        ExtractChrinfo("input/input.sorted.bam"), "output/")
```

The outputs are a combined annotated APA clusters file and three files (genes.tsv, barcodes.tsv and martrix.mtx) in a directory named "matrix/" under "output/". 

The first 8 columns in combined annotated APA clusters file are the same as the output of *AnnoAPAc()*, and the last four columns are cluster unique id, cluster id in each gene, cluster count per gene and min-max normalization value of location in each gene.

The matrix can be easily loaded into [Seurat](https://satijalab.org/seurat/):

```R
library(Seurat)
pbmc.data  <- Read10X("output/matrix/")
```



## Sample 🎥

The [3k PBMC dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.0.0/pbmc3k) from 10x Genomics were used for testing.

The sample outputs can be found [here](https://github.com/HHengUG/hatpal_example/tree/main).

A script to run hatpal automatically can be found [here](https://github.com/HHengUG/hatpal_example) .



## Contact 📨

Question and feedback are welcome at:

HHeng: hengheng at genomics dot cn



## News ✨

2021/05/03

- update new version (https://github.com/HHengUG/hatpal2)

2021/01/20

- update README.md

2020/11/23

- update README.md
- update samples

2020/11/6

- update README.md

2020/11/5: New releases

- hatpal 0.9.1

