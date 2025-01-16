# 甲基化流程
ACE-seq方法原理：   
A3A酶将未修饰的胞嘧啶C和5mC转换为尿嘧啶U，5hmC保持不变，并使用β-葡萄糖基转移酶（βGT）将5hmC转化为5-葡萄糖基羟甲基胞嘧啶（5ghmC）以保护5hmC。    
重亚硫酸盐测序BS-seq原理：   
亚硫酸盐将未修饰的C转换为U，PCR后变成T，5mC和5hmC不变。
- 工具：    
aria2:1.37.0  
fastqc:0.12.1  
trimgalore:0.6.10  
bismark:0.24.2   
## 数据准备：
    测序数据GEO号：[GSE116016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116016)   
    使用其中部分测序数据：野生型SRR7368841、SRR7368842；TetTKO SRR7368845    
    小鼠基因组：[Ensembl](https://www.ensembl.org/info/data/ftp/index.html)    

- 测序数据下载  
ENSEMBL查询GEO号后`Run Selector`中SRR号，再到ENA中用SRR号获取ftp地址，用`aria2`下载
```
cd ~/project/mouse/sequence
aria2c -d . -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/001/SRR7368841/SRR7368841.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/002/SRR7368842/SRR7368842.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/005/SRR7368845/SRR7368845.fastq.gz
```
ps:从NCBI用sratoolkit中`prefetch`下载后的`.sra`文件后，用`md5sum`命令计算md5码，可以与ENA页面的`fastq_md5`码对比，确认数据是否完整
- 参考基因组数据、
ENSEMBL下载
```
cd ~/project/mouse/genome
wget ftp://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz
gzip -d Mus_musculus.GRCm39.dna.toplevel.fa.gz
mv Mus_musculus.GRCm39.dna.toplevel.fa.gz GRCm39.fa
```
## 质控清洗
使用`fastqc`和`trimgalore`：
```
cd ~/project/mouse/sequence
fastqc -t 3 -o ../output/fastqc *.gz
trim_galore -o ../output/trim --fastqc *.fastq.gz
```
`--fastqc`参数：trim完直接再fastqc质控一次

## 甲基化分析
使用`Bismark`工具，该工具基于BE-seq，也可以用于ACE-seq，以19号染色体为例
### 基因组索引(亚硫酸盐版本)
使用`Bowtie`和`Bowtie2`建立索引
```
bismark_genome_preparation --bowtie2 ~/project/mouse/genome_chr19
```
### **比对**
甲基化分析核心步骤
```
mkdir ~/project/mouse/output/bismark_align_chr19
cd ~/project/mouse/sequence

# 比对
bismark -o ../output/bismark_align_chr19 --parallel 4 --genome_folder ../genome_chr19/ *.fastq.gz
# 合并WT两个文件
cd ~/project/mouse/output/bismark_align_chr19
samtools cat -o SRX4241790_trimmed_bismark_bt2.bam SRR7368841_bismark_bt2.bam SRR7368842_bismark_bt2.bam
```

### 比对reads去除重复
使用`deduplicate_bismark`
```
cd ~/project/mouse/output
mkdir deduplicate
cd ~/project/mouse/output/bismark_align_chr19

deduplicate_bismark --bam --output_dir ../deduplicate SRR7368845_bismark_bt2.bam SRX4241790_trimmed_bismark_bt2.bam
```

### 提取甲基化信息
使用`bismark_methylation_extractor`从比对结果提取甲基化信息
```
cd ~/project/mouse/output
mkdir methylation_extractor
bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \
    --cytosine_report --genome_folder ../genome_chr19 \
    -o ./methylation_extractor ./deduplicate/*.bam
```
`--bedGraph`：提取完甲基化信息后，甲基化输出写入一个排序后的bedGraph文件中，该文件包含胞嘧啶位置及甲基化状态。
`--cytosine_report`：bedGraph文件转换后，使用该选项生成全基因组范围内所有胞嘧啶的甲基化报告。

## 下游分析
下游分析，包括寻找特定位点、检测差异化甲级位点DML、差异化甲基区域DMR，使用R包`DSS`差异甲基化分析。
### 输入数据准备
DSS要求每个CG位点上总结为以下信息：染色体编号、基因组坐标、总reads数、甲基化的reads数   
所需输入数据从bismark结果中的`.cov`文件转换，提取`.cov`文件中的列，生成DSS要求的输入格式   
`.cov`文件包含的列：`chr`、`start`、`end`、`methylation`、`count methylated`、`count unmethylated`   
```
cd ~/project/mouse/output
mkdir Ranalysis
cd Ranalysis

# 存储数据路径
cov_WT_path="$HOME/project/mouse/output/methylation_extractor/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
cov_TetTKO_path="$HOME/project/mouse/output/methylation_extractor/SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz"
# 复制到Ranalysis文件
cp $cov_WT_path .
cp $cov_TetTKO_path .
# 解压
gzip -d SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov.gz
# 转换.cov格式为.txt格式
cp SRR7368845_bismark_bt2.deduplicated.bismark.cov SRR7368845_bismark_bt2.deduplicated.bismark.txt
cp SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.txt
```
```R
library(tidyr)
library(dplyr)

file.names <- c("./SRX4241790_methylation_result.txt", "./SRR7368845_methylation_result.txt")

func_read_file <- function(file_name){
	dir_vec <- strsplit(file_name, split = "/")[[1]]
	len <- length(dir_vec)
	file_prefix = substring(dir_vec[len], 0, nchar(dir_vec[len]) - 4)
	file_save_path = substring(file_name, 0, nchar(file_name) - nchar(dir_vec[len]))
	print(paste("File", file_name, "is being importing and this may take a while..."), sep = "")
	rawdata_df <- read.table(file_name, header = F, stringsAsFactors = F)
	print("Importing file is finished!")
	colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
	write.table(rawdata_df, paste(file_save_path, file_prefix, "_transfered.txt", sep = ""), row.names = F )
}

lapply(file_names, func_read_file)
```
