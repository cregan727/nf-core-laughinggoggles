#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*

Cellsnp-lite modules for all varients of running cellsnp-lite

*/


// cellsnp-lite for 10x - used in all subworkflows

process cellsnp_lite_10x {

    /*

	This process takes as input a single 10x bam file
	sample name, barcodes list, and a regionvcf and does
	variant calling ahead of reference based or reference
	free demultiplexing with vireo

    */



  label 'process_high'

  publishDir params.publishDir, mode: 'copy'

  cpus 20

container {
    container = 'jeffverboon/cellsnplite:latest'
}

  input:
    tuple val(sample),
	file(bam),
        file(bai),
	file(barcodes),
    file(regionvcf)

  output:
    path "${sample}_outdir/"

  script:
    """
cellsnp-lite \
        -s $bam \
        -b $barcodes \
        -O ${sample}_outdir/ \
        -R $params.regionvcf \
        -p 20 \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip
    """


}


// cellsnp-lite bamfile input


process cellsnp_lite_bams {

    /*

        This process takes as input a list of STAR mapped bams             
        from an RNA-seq experiment, a list of sample names, and
	a regionvcf and does variant calling ahead of reference
	based or reference free demultiplexing with vireo

    */



  label 'process_high'

  cpus 20

container {
    container = 'jeffverboon/cellsnplite:latest'
}

  input:
    path bamdir
    file params.regionvcf

  output:
    file "outdir/cellSNP.cells.vcf.gz"

  script:
    """
ls bamdir/*

cellsnp-lite \
        -s \$(head -n 1 bamdir/bamsamplesheet.txt) \
        -I \$(tail -n 1 bamdir/bamsamplesheet.txt) \
        -O outdir/ \
        -R $params.regionvcf \
        -p 20 \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --cellTAG None \
        --UMItag None \
        --genotype \
        --gzip
    """

}



// cellsnp-lite jplate input

