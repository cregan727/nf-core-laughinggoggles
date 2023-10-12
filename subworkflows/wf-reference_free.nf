include { cellsnp_lite_10x } from '../modules/cellsnplite.nf'
include { vireo_reffree } from '../modules/vireo.nf'


input_data = Channel.fromPath(params.samplesheet)
                      .splitCsv(header: true, sep: ',')


workflow reffree_workflow {

    /*

	This version of the workflow takes a samplesheet defining the 10x samples
	and an expected number of samples.

    */

  // take input samplesheet, for each row perform cellsnp-lite and then reffree vireo
input_data.map { row ->
    tuple(row['Sample'], file(row['bam']), file(row['bam'] + ".bai"), file(row['barcode']), file(params.regionvcf))
} | cellsnp_lite_10x | vireo_reffree | view

}
