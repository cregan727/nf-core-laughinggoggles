# nf-decoder-ring

Internal WIP nextflow pipeline for Genetic Demultiplexing of 10x Genomics libraries

## Example Usage (reference free version):

```
nextflow run cregan727/nf-decoder-ring -r main --samplesheet $PWD/samples.analysis.csv \
        --regionvcf "path/to/vcf.gz"  \
        --publishDir $PWD/published_results \
        --workflow "ref_free" --numsamples 2 \
        -profile singularity -resume 
```

When the samples.analysis.csv looks like:

```
Sample,bam,barcode
Sample1,/path/to/Projectdir/count/Sample1/outs/possorted_genome_bam.bam,/path/to/Projectdir/count/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
Sample2,/path/to/Projectdir/count/Sample2/outs/possorted_genome_bam.bam,/path/to/Projectdir/count/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
```

__________________________________________________________________


## Example Usage (reference bulk RNA-seq version):

```
nextflow run cregan727/nf-decoder-ring -r main --samplesheet $PWD/samples.analysis.csv \
        --regionvcf "path/to/vcf.gz"  \
        --publishDir $PWD/published_results \
        --workflow "ref_bams" \
        --samplesheet_bams $PWD/bams_samplesheet.csv \
        -profile singularity -resume
```

When the samples.analysis.csv is the same and the bams_samplesheet.csv looks like (and the directory containing the bam files also contains a .bam.bai file):

```
/path/to/Aligned.sortedByCoord.out.bam,Sample1
/path/to/other/Aligned.sortedByCoord.out.bam,Sample2
```
__________________________________________________________________

## Example Usage (reference plate based bulk RNA-seq):

```
nextflow run cregan727/nf-decoder-ring -r main --samplesheet $PWD/samples.analysis.csv \
        --regionvcf "path/to/vcf.gz"   \
        --publishDir $PWD/published_results \
        --workflow "ref_plate" \
        --samplesheet_plate $PWD/jplate_samplesheet.csv \
        -profile singularity -resume 
```

When the samples.analysis.csv is the same and the jplate_samplesheet.csv looks like (and the directory containing the bam files also contains a .bam.bai file):

```
/path/to/Aligned.sortedByCoord.out.bam,/path/to/barcode.csv
```
And the barcode csv can either be a 10x style barcode.tsv.gz file or like this:
```
ATCGTAGC,Sample1
TCGGTCTA,Sample2
```


____________________________________________________________________
## A Note on Sequencing Depth

While we have not done extensive testing on how deeply to sequence your reference samples or single cell libraries for ideal genetic demutiplexing, if you have a large number of unassigned samples and a lower than expected number of doublets it may be a sign that you should sequence your reference and/or single cell library deeper to increase confidence in the vireo calls. Here is an example of an experiment in which the reference libraries were sequenced to an average of 5M reads per cell and the 10x 3' V2 libraries were sequenced to 56,000 reads per cell and downsampled ahead of cellsnp-vireo with the pipeline described above. When only a small number of reference reads were used (< 1 M reads per sample) you can see the number of unassigned samples is much larger and the number of doublets is pretty low, but as the sequencing depth increases more cells get assigned to their samples and more doublets get identified. Similarly once the number of reads per cell for the single cell library dropped below 10,000 the number of unassigned cells increased dramatically while the number of doublets decreased.

![How_deep_should_I_seq_plates](https://github.com/cregan727/nf-laughing-goggles/assets/68451521/8014d16b-3ba5-48c2-84d4-df398c6df2f7)
![how_deep_to_seq_scrna](https://github.com/cregan727/nf-laughing-goggles/assets/68451521/688514f6-81a2-4580-91ec-5c3659fe73cc)


## Other  

Where can I get a set of SNPs to use for CellSNP? Check out this resource provided in the original python version of CellSNP: https://sourceforge.net/projects/cellsnp/files/SNPlist/ 

### Acknowledgments  

Xianjie Huang, Yuanhua Huang, Cellsnp-lite: an efficient tool for genotyping single cells, Bioinformatics, 2021;, btab358, https://doi.org/10.1093/bioinformatics/btab358  
Yuanhua Huang, Davis J. McCarthy, and Oliver Stegle. Vireo: Bayesian demultiplexing of pooled single-cell RNA-seq data without genotype reference. Genome Biology 20, 273 (2019)
