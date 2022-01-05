# Fec
Fec is a error correction tool based on two-rounds overlapping and cacheing. The first round overlapping will find a number of overlaps quickly. Based on the overlaps, some reads can be corrected immediately, and the rest reads will be performed the second round overlapping sensititively to find as more overlaps as possible.
## Dependencies
A compiler that supports C++11 is needed to build nanopolish. Development of the code is performed using [gcc-10.1](https://gcc.gnu.org/gcc-10/).
- [minimap2-c9874e2](https://github.com/lh3/minimap2/)
## Installation instructions
- ### Install Fec from github
You can download and compile the latest code from github as follows:
```
git clone https://github.com/zhangjuncsu/Fec.git
cd Fec
make
cd
```
After installation, all the executables can be found in Fec/Linux-amd64/bin. The folder name Linux-amd64 will vary in operating system.
- ### Install minimap2 from github
```
git clone https://github.com/lh3/minimap2.git
cd minimap2
git checkout c9874e2
make
export PATH=pwd:$PATH
```
## Usage
For small datasets:
```
minimap2 -x ava-pb -w 20 -K 2g -t $THREADS reads.fq reads.fq | awk '{ if($4 - $3 >= 0.5 * $2 || $9 - $8 >= 0.5 * $7) print $0}' > ovlp.paf
/path/to/Fec/Linux-amd64/bin/fec -t $THREADS -r 0.6 -a 1000 -c 4 -l 2000 ovlp.paf reads.fq corrected.fasta
```
For large datasets such as human:
```
minimap2 -x ava-pb -w 20 -K 2g -f 0.005 -t $THREADS reads.fq reads.fq | awk '{ if($4 - $3 >= 0.2 * $2 || $9 - $8 >= 0.2 * $7) print $0}' > ovlp.paf
/path/to/Fec/Linux-amd64/bin/fec -t $THREADS -r 0.6 -a 1000 -c 4 -l 2000 -m 0.005 -f 0.2 ovlp.paf reads.fq corrected.fasta
```
Fec now only works with uncompressed FASTA or FASTQ formats for reads. [Minimap2](https://github.com/lh3/minimap2) is used to perform overlapping. As minimap2 will find some short overlaps, Fec use **awk** to filter those overlaps. For large datasets, such as huamn daasets, you can set option **-f 0.005** to speedup overlapping. The second round overlapping will be performed internally in correcting process, you can use option **-m** to filter out top FLOAT fraction of repetitive minimizers. and you can use option **-f** to filter out short overlaps found in the second round overlapping.