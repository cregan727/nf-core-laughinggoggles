# nf-core-laughinggoggles

Internal WIP nextflow pipeline for Genetic Demultiplexing of 10x Genomics libraries

## Example Usage (reference free version):

```
nextflow run cregan727/nf-core-laughinggoggles -r main --samplesheet $PWD/samples.analysis.csv \
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
nextflow run cregan727/nf-core-laughinggoggles -r main --samplesheet $PWD/samples.analysis.csv \
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
nextflow run cregan727/nf-core-laughinggoggles -r main --samplesheet $PWD/samples.analysis.csv \
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
