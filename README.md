# 甲基化流程
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
使用`Bismark`工具，该工具基于BE-seq，也可以用ACE-seq
### 基因组索引(亚硫酸盐版本)
使用`Bowtie`和`Bowtie2`建立索引
```
bismark_genome_preparation --bowtie2 ~/project/mouse/genome
```
### **比对**
甲基化分析核心步骤
```
mkdir ~/project/mouse/output/bismark_align
cd ~/project/mouse/sequence

# 比对
bismark -o ../output/bismark_align --parallel 4 --genome_folder ../genome *.fastq.gz
# 合并
samtools cat -o SRX4241790