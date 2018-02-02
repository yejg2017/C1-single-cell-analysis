# C1 single cell analysis项目（帮别人做的作业，哈哈，just a homework，haha）
该项目主要分析**Human,Monkey,Mouse**三个数据集的Genecount数据集

## 1）数据说明：
这三个数据集都是genecounts数据，G x N的dataframe格式（当然是用R处理整合了raw data），G表示基因，N表示样本.
样本的名字格式：**mkc005_6um_H3_S275_R1_001**都是类似这样，譬如，*mkc005表示个体来源，6um表示细胞的大小*

## 2）数据处理：
### 2.1) Human:
* 删除hc001,shoutiao和细胞大小除6um,10um,20um的样本;
* 计算剩余每个样本的全部基因counts的总和，得到所有样本的all genes total counts,然后对total counts做16分位数，去掉1/16,15/16两端的极端,根据这标准删除样本,此时大概一共删除了50个样本;
* 在构建以下Seurat object时，取标准min.cells,min.genes分别为10,2,即分别代表每个细胞至少包含10个表达基因,每个基因至少在2个基因中表达.
* 完成以上处理步骤，数据从56505个基因,1099个样本变为了18880个基因和895个样本.

### 2.2) Mouse,Monkey:
* 计算剩余每个样本的全部基因counts的总和，得到所有样本的all genes total counts,然后对total counts做16分位数，去掉1/16,15/16两端的极端,根据这标准删除样本.
* 在构建以下Seurat object时，取标准min.cells,min.genes分别为10,2,即分别代表每个细胞至少包含10个表达基因,每个基因至少在2个基因中表达.

以上数据处理都是根据自己的一个经验测试的而已，没有很大的根据，如有更好的方法，[please contact me email]:1772898816@qq.com or yejg2013@gmail.com


## 3)使用的工具 **rstudio**
主要使用的 r packages：
``` r
library(Seurat)
library(data.table) 
library(NMF)//r
library(rsvd)//r
library(Rtsne)//r
library(ggplot2)
library(cowplot)
library(sva)
library(igraph)
library(cccd)
library(KernSmooth)
library(beeswarm)
library(stringr)
library(formatR)
source('../tools.R')  #自己写的一些工具函数
library(DESeq2)
```
不一一介绍安装（自己google，度娘去）,**Seurat,ggplot2**是主要的分析工具，重要的事情说三遍

## 4)结果展示
结果都以*markdown*的格式展示在各自文件夹**Human.Monkey,Mouse**的*docx,html*文件中

## 5）示例图片







