include { cellsnp_lite_10x } from '../modules/cellsnplite.nf'
include { cellsnp_lite_plate } from '../modules/cellsnplite.nf'
include { vireo_ref } from '../modules/vireo.nf'


// load in Single Cell Data

input_data = Channel.fromPath(params.samplesheet)
                      .splitCsv(header: true, sep: ',')

workflow ref_plate {
    /*

        This version of the workflow takes a samplesheet defining the 10x samples 
        and another samplesheet defining the path to the bams and the names of 
        those samples for use as a reference to genetically demultiplex. 

    */



    // take input as bam files and run cellsnp-lite
    //preprocess_output = process_bam_samplesheet(params.samplesheet_bams)
    //preprocess_output.view()
    input_ref = Channel.fromPath(params.samplesheet_plate)
                            .splitCsv(header: false, sep: ',')
    bam_cellsnp = input_ref.map { sheet ->
    tuple(file(sheet[0]), file(sheet[0] + ".bai"), file(sheet[1]), file(params.regionvcf))
} | cellsnp_lite_plate


    //perform cellsnp-lite on the 10x bams provided in the samplesheet
    wf1_out = input_data.map { row ->
    tuple(row['Sample'], file(row['bam']), file(row['bam'] + ".bai"), file(row['barcode']), file(params.regionvcf))
} | cellsnp_lite_10x



    //run reference based vireo
    //wf1_out.each {path -> vireo_ref(bam_cellsnp, path)} 
    scatter (path in wf1_out) {
        vireo_ref(bam_cellsnp, path)
    }

}
