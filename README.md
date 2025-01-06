# 甲基化流程
### 数据准备：
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
```