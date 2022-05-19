# Fec
Fec is a error correction tool based on two-rounds overlapping and caching. The first round overlapping will find a number of overlaps quickly. Based on the overlaps, some reads can be corrected immediately, and the rest reads will be performed the second-round overlapping using finely tuned to find as more overlaps as possible.
## Dependencies
A compiler that supports C++11 is needed to build Fec. Development of the code is performed using [gcc-10.1](https://gcc.gnu.org/gcc-10/).
- [minimap2](https://github.com/lh3/minimap2/)
## Installation instructions
- ### Install Fec from bioconda
You can install Fec from bioconda channel
```
conda install -c bioconda fec
```
- ### Install Fec from github
You can download and compile the latest code from github as follows:
```
git clone https://github.com/zhangjuncsu/Fec.git
cd Fec
make
export PATH=`pwd`/Linux-amd64/bin:$PATH
cd
```
After installation, all the executables can be found in ./Linux-amd64/bin. The folder name Linux-amd64 will vary in operating system.
- ### Install minimap2 from github
```
git clone https://github.com/lh3/minimap2.git
cd minimap2
make
export PATH=`pwd`:$PATH
```
## Usage
**For small PacBio datasets:**
```
minimap2 -x ava-pb -w 20 -K 2g -t $THREADS reads.fq reads.fq | awk '{ if($4 - $3 >= 0.5 * $2 || $9 - $8 >= 0.5 * $7) print $0}' > ovlp.paf
Fec -t $THREADS -r 0.6 -a 1000 -c 4 -l 2000 ovlp.paf reads.fq corrected.fasta
```
**For large PacBio datasets such as human:**
```
minimap2 -x ava-pb -w 20 -K 2g -f 0.005 -t $THREADS reads.fq reads.fq | awk '{ if($4 - $3 >= 0.2 * $2 || $9 - $8 >= 0.2 * $7) print $0}' > ovlp.paf
Fec -t $THREADS -r 0.6 -a 1000 -c 4 -l 2000 -m 0.005 -f 0.2 ovlp.paf reads.fq corrected.fasta
```

**For nanopore datasets:**
We use two round correction for better accuracy.
```
minimap2 -x ava-ont -w 20 -K 2g -f 0.005 -t $THREADS reads.fq reads.fq | awk '{ if($4 - $3 >= 0.2 * $2 || $9 - $8 >= 0.2 * $7) print $0}' > ovlp.paf
Fec -x 1 -t $THREADS -r 0.6 -a 400 -c 0 -l 1000 -m 0.005 -f 0.2 ovlp.paf reads.fq corrected1.fasta
minimap2 -x ava-ont -w 20 -K 2g -f 0.005 -t $THREADS corrected1.fasta corrected1.fasta | awk '{ if($4 - $3 >= 0.2 * $2 || $9 - $8 >= 0.2 * $7) print $0}' > ovlp.paf
Fec -x 1 -R -t $THREADS -r 0.6 -a 1000 -c 4 -l 2000 -m 0.005 -f 0.2 ovlp.paf corrected1.fasta corrected2.fasta
```
Fec now only works with uncompressed FASTA or FASTQ formats for reads. [Minimap2](https://github.com/lh3/minimap2) is used to perform overlapping. As minimap2 will find some short overlaps, Fec use **awk** to filter those overlaps. For large datasets, such as human daasets, you can set minimap2 option **-f 0.005** to speedup overlapping. The second round overlapping will be performed internally in correcting process, you can use option **-m** to filter out top FLOAT fraction of repetitive minimizers, and use option **-f** to filter out short overlaps found in the second round overlapping.