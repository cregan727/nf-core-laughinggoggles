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
    container = 'ghcr.io/cregan727/decoder_ring_cellsnplite:latest'
}

  input:
    tuple val(sample),
	file(bam),
        file(bai),
	file(barcodes),
    file(regionvcf)
    params.minUMI

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
        --minCOUNT $params.minUMI \
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
    container = 'ghcr.io/cregan727/decoder_ring_cellsnplite:latest'
}

  input:
    path bamdir
    file params.regionvcf
    params.minUMI

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
        --minCOUNT $params.minUMI \
        --cellTAG None \
        --UMItag None \
        --genotype \
        --gzip
    """

}



// cellsnp-lite jplate input


process cellsnp_lite_plate {

    /*

        This process takes as input a single STAR mapped bam
        from a plate, barcodes list, and a regionvcf and does
        variant calling ahead of reference based demultiplexing
        with vireo

    */



  label 'process_high'

  publishDir params.publishDir, mode: 'copy'

  cpus 20

container {
    container = 'ghcr.io/cregan727/decoder_ring_cellsnplite:latest'
}

  input:
    tuple file(bam),
        file(bai),
        file(barcodes),
    file(regionvcf)
    params.minUMI

  output:
    file "outdir/cellSNP.cells.vcf.gz"

  script:
    """
cut -d ',' -f 1 $barcodes > barcodes_for_cellsnp.csv

cellsnp-lite \
	-s $bam \
        -b barcodes_for_cellsnp.csv \
        -O outdir/ \
        -R $params.regionvcf \
        -p 20 \
        --minMAF 0.1 \
        --minCOUNT $params.minUMI \
        --genotype \
        --gzip
        
declare -A barcode_to_sample

while IFS=',' read -r barcode sample; do
  if [[ ! -z "\$sample" ]]; then
    barcode_to_sample["\$barcode"]="\$sample"
  else
    echo "No sample names found. Returning cellSNP.cells.vcf.gz with barcodes"
    echo "If this was not what you expected confirm there are no blank lines in your barcodes file"
    echo "and there is a sample name for every barcode in your barcodes file has a sample name"
    exit 0  # Exit with a non-error exit code
  fi
done < $barcodes

  for barcode in "\${!barcode_to_sample[@]}"; do
    sample="\${barcode_to_sample[\$barcode]}"
    gunzip -c outdir/cellSNP.cells.vcf.gz | sed "s/\$barcode/\$sample/g" | gzip > outdir/cellSNP.cells_modified.vcf.gz
    mv outdir/cellSNP.cells_modified.vcf.gz outdir/cellSNP.cells.vcf.gz
  done


    """


}

