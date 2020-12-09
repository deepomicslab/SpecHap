# SpecHap

SpecHap is an ultra fast phasing algorithm based on spectral graph analysis. SpecHap currently support general WGS sequencing, Hi-C and 10X-linked reads.


## Getting Started

To build SpecHap, run the following command:

```
cd /path/to/SpecHap/src
mkdir build
cd build
cmake ..
make && make install
```

To see the option SpecHap support, run the following command:

```
SpecHap --help
```

A modified utility software ExtractHair, originally from [HAPCUT2](https://github.com/vibansal/HapCUT2), is needed for data preprocessing. To install, run
```
cd /path/to/SpecHap/
git submodule init && git submodule update
cd hair-src
make && make install
```
### Prerequisites

SpecHap relies on ARPACK for Eigen-calculation. To gain stable utilization, we recommend [arapack-ng](https://github.com/opencollab/arpack-ng).

Arpack-ng can be easily compile with cmake, ensure you have BLAS and LAPACK installed before compiling.

To build [arapack-ng](https://github.com/opencollab/arpack-ng), try run
```
cd /path/to/arpack-ng/
sh bootstrap
./configure --enable-icb
make
make && make install
```

[Htslib](https://github.com/samtools/htslib) is also required. To install [htslib](https://github.com/samtools/htslib), simply run
```
cd /path/to/htslib
autoheader  #required if htslib is cloned from github
autoconf    #required if htslib is cloned from github
./configure
make && make install
```



### Using SpecHap

#### Data preprocessing
SpecHap requires at least a fragment file and a bgziped and indexed VCF to perform phasing. Ensure your VCF is sorted by position. The fragment file is also required to be sorted, the SpecHap is set to sort the fragment for you automatically on default option.

To generate the fragment file, run the following command
```
extractHAIRS --bam /your/bam/file --VCF /your/vcf/file --out fragment_file
``` 
With Hi-C sequenced file, try 
```
extractHAIRS --bam /your/bam/file --VCF /your/vcf/file --out fragment_file --hic 1
```
With 10X linked reads, try 
```
extractHAIRS --bam /your/bam/file --VCF /your/vcf/file --out fragment_file --10x 1
```
You also need a bed file indicating each barcode's inferred spanning range. You can use the BarcodeExtract to do your job
```asm
BarcodeExtract /you/bam/file barcode_spnanning.bed
bgzip -c barcode_spanning.bed > barcode_spanning.bed.gz
tabix -p bed barcode_spanning.bed.gz
```

#### Run SpecHap
The detailed usage can be found by run
```asm
SpecHap --help
```

Notice that for 10X linked reads, we recommend you to run with window size 200 and overlap length 60

### Author
SpecHap is developed by Delta team under the supervision of Dr. Li Shuaicheng, City University of Hong Kong

To contact us, send email to [yonghanyu2@cityu.edu.hk](yonghanyu2@cityu.edu.hk)

## Built With

* [htslib](https://github.com/samtools/htslib)
* [arapack-ng](https://github.com/opencollab/arpack-ng)
* [Eigen3](https://eigen.tuxfamily.org/dox/)
* [Lean Mean C++ Option Parser](https://eigen.tuxfamily.org/dox/)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


