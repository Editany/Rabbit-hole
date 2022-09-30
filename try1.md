## $\color{#8C0044}{寻找基因突变位点-基本流程}$

### 一些基本操作

##### 查看服务器情况 
top
htop
按q退出

##### 用screen防止断连，screen把程序给了后台，不会中断。exit掉screen程序也不会断（只是不能再控制）
screen -S xxx 创建
screen -ls 查看已有screen
screen -r 打开某个screen
按ctrl+A+D关闭screen
exit或按ctrl+D删除screen

##### bash脚本
touch bash.sh 创建
chmod u+x bash.sh 改权限
vim bash.sh 开始编辑
按I开始编辑（第一次进入不用按，会直接进入INSERT模式）
按Esc，之后输入：wq保存退出
bash bash.sh运行脚本


##### -----------------------------------------------------------------------------------------------

* 以下过程使用了两个表型不同的玉米蛇样本（PAGU009、PAGU044），实际研究中寻找性状对应的基因位点会需要20-30个左右的样本。

##### -----------------------------------------------------------------------------------------------

### 一、bwa比对，生成sam文件

##### 1.建立bwa索引
fasta文件从ncbi里下载（搜索物种名，下载基因组fasta），下载后上传到服务器中。
![98afd4f1c7e1f2120b9c140fde3581cd.png](en-resource://database/510:1)

bwa index 文件在服务器中的路径 建立索引
```
bwa index PanGut3.0/GCF_001185365.1/GCF_001185365.1_UNIGE_PanGut_3.0_genomic.fna
```
#这一步用时较短，玉米蛇基因组用时约2h


##### 2.bwa比对

```
bwa mem [options] ref.fa read1.fq read2.fq > mem-paired-ends.sam
```
这一步将生成sam格式的比对文件
详细内容可以戳这里[sam文件格式介绍](http://samtools.github.io/hts-specs/SAMv1.pdf)
sam格式下一段序列的比对信息有两部分，header section和alignment section。
头部区包括@HD格式版本和排序类型、@SQ参考序列信息和长度、@RG样品信息、@PG比对软件及版本

代码中需要使用到的options如下
-t 线程数
-M 将shorter split hits标记为次优
-R 设置reads标头，放到一对单引号中，会成为sam文件中的RG部分（序列说明），形如：@RG ID:xxx SM:xxx LB:xxx PU:xxx PL:xxx。写入时用\t制表符分割，例如这里引号内的内容为 '@RG\tID:PAGU009\tPL:ILLUMINA\tLB:PAGU009\tSM:PAGU009'


ID：Read grope的ID号,每个RG都有一个独特的ID
PL：测序平台，可选CAPILLARY,DNBSEQ(MGI/BGI),ELEMENT,HELICOS,ILLUMINA,IONTORRENT,LS454,ONT,PACBIO,SOLID,ULTIMA
LB：文库名，文库是许多DNA片段两端接上特定接头形成的DNA混合物
SM：样品名，表示该read group所存储的reads所来源的样本名称
PU不必须：测序平台单元，例如Illumina是{flowcell}.{lane}.{barcode}
flowcell流动槽 lane泳道 barcode混合样品中每个样品的身份证

ref.fa是已经建立了索引的参考基因组
read1.fq read2.fq是样本测序结果文件夹里的两个文件
">"后是设定的输出文件名xxx.sam


```
bwa mem -t 5 -M -R '@RG\tID:PAGU009\tPL:ILLUMINA\tLB:PAGU009\tSM:PAGU009' PanGut3.0/GCF_001185365.1/GCF_001185365.1_UNIGE_PanGut_3.0_genomic.fna snake_sample/PAGU/PAGU-009_FDSW192264952-1a/PAGU-009_FDSW192264952-1a_1.clean.fq.gz snake_sample/PAGU/PAGU-009_FDSW192264952-1a/PAGU-009_FDSW192264952-1a_2.clean.fq.gz > PAGU009.sam

bwa mem -t 5 -M -R '@RG\tID:PAGU044\tPL:ILLUMINA\tLB:PAGU044\tSM:PAGU044' PanGut3.0/GCF_001185365.1/GCF_001185365.1_UNIGE_PanGut_3.0_genomic.fna snake_sample/PAGU/PAGU-044_FDSW192264974-1a/PAGU-044_FDSW192264974-1a_1.clean.fq.gz snake_sample/PAGU/PAGU-044_FDSW192264974-1a/PAGU-044_FDSW192264974-1a_2.clean.fq.gz > PAGU044.sam
```
#生成约70G的sam文件用时约12h

##### -----------------------------------------------------------------------------------------------

### 二、sam转bam及根据FLAG、MAPQ过滤
#这几步用时都比较短

##### 1.samtools view 转换sam与bam文件，并对bam文件进行过滤操作
bam文件是sam文件的二进制格式，占据内存较小且运算速度快
```

samtools view [options] in.sam
```
使用的options如下
-h 输出文件中带header
-f 保留该FLAG值的序列 
-F 保留除该FLAG值的序列
-q 只输出MAPQ值大于该值的序列
-b 输出格式为bam（默认输出为sam）
-S 输入格式为sam（默认输入为bam） 
-o 输出文件名

sam及bam文件中的比对区第二列为FLAG信息，第五列为MAPQ（mapping quality）

FLAG表现为多个数字的加和，这里3=1+2，3852=4+8+256+512+1024+2048
不同数字代表不同的含义
1  该read是成对的paired reads中的一个 \
2  paired reads中每个都正确比对到参考序列上 \
4  该read没比对到参考序列上 \
8  与该read成对的matepair read没有比对到参考序列上 \
16  该read其反向互补序列能够比对到参考序列 \
32  与该read成对的matepair read其反向互补序列能够比对到参考序列 \
64  在paired reads中，该read是与参考序列比对的第一条 \
128 在paired reads中，该read是与参考序列比对的第二条 \
256  该read是次优的比对结果 \
512  read没有通过质量控制 \
1024  由于PCR或测序错误产生的重复reads \
2048  补充匹配的read

MAPQ 即mapping的错误率的-10log10值，0-60，0为unmapped，一般到60代表比对几乎完全了，255表示值不可用
```
samtools view -h -f 3 -F 3852 -q 30 -b -o PAGU009.bam PAGU009.sam
samtools view -h -f 3 -F 3852 -q 30 -b -o PAGU044.bam PAGU044.sam
```

* tip 一个查看进度的方法：查看文件行数

```
samtools view -c xxx.bam
```


##### 2.去除重复序列

① samtools sort排序
```
samtools sort [options] [-o out.bam] [in.sam|in.bam|in.cram]
```
使用的options如下
-n 按照【测序序列名称】（QNAME，sam格式比对部分第一列）排序
-o 设置输出文件名称
```
samtools sort -n -o PAGU009.namesort.bam PAGU009.bam
samtools sort -n -o PAGU044.namesort.bam PAGU044.bam
```

 

② fixmate为【以名称排序】的bam在最后的列填入配对坐标，ISIZE（推测的插入序列大小）和配对相应的flag
```
samtools fixmate [-rpcmu] [-O format] in.namesort.bam out.bam
```
使用的options如下
-m 添加mate score tags，后续markdup选择最佳read时需要
```
samtools fixmate -m PAGU009.namesort.bam PAGU009.fixmate.bam
samtools fixmate -m PAGU044.namesort.bam PAGU044.fixmate.bam
```

③ 再次排序（samtools sort）
 不加-n 以【最左端比对位置坐标】（POS，sam格式比对部分第四列）排序
-m XX 每个线程最多消耗XX内存，阻止创建过多tmp文件
```
samtools sort -m 4G -o PAGU009.positionSort.bam PAGU009.fixmate.bam
samtools sort -m 4G -o PAGU044positionSort.bam PAGU004.fixmate.bam
```

④ markdup 在【以坐标排序】的文件中标记重复的比对
-r 去除重复reads
```
samtools markdup -r PAGU009.positionSort.bam PAGU009.markdup.bam
samtools markdup -r PAGU044.positionSort.bam PAGU044.markdup.bam
```

##### -----------------------------------------------------------------------------------------------
