include { cellsnp_lite_bams } from '../modules/cellsnplite.nf'
include { cellsnp_lite_10x } from '../modules/cellsnplite.nf'
include { vireo_ref } from '../modules/vireo.nf'
include { process_bam_samplesheet } from '../modules/other.nf'


input_data = Channel.fromPath(params.samplesheet)
                      .splitCsv(header: true, sep: ',')


workflow ref_bams {
    /*

	This version of the workflow takes a samplesheet defining the 10x samples 
	and another samplesheet defining the path to the bams and the names of 
	those samples for use as a reference to genetically demultiplex. 

    */


    // take input as bam files and run cellsnp-lite
    preprocess_output = process_bam_samplesheet(params.samplesheet_bams)
    bam_cellsnp = cellsnp_lite_bams(preprocess_output, file(params.regionvcf))


    //perform cellsnp-lite on the 10x bams provided in the samplesheet
    wf1_out = input_data.map { row ->
    tuple(row['Sample'], file(row['bam']), file(row['bam'] + ".bai"), file(row['barcode']), file(params.regionvcf))
} | cellsnp_lite_10x

    //run reference based vireo
    wf1_out.each {path -> vireo_ref(bam_cellsnp, path)} 


}

