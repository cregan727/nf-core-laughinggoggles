/*

All of the vireo options - either reference based or reference free

*/



// reference free vireo

process vireo_reffree {

    label 'process_high'
       
    publishDir params.publishDir, mode: 'copy'
    
    container {
    container = 'ghcr.io/cregan727/decoder_ring_vireo:latest'
    }
    
    input:
    path sample_outdir
    params.numsamples

    output:
    path "*/vireo_results"
    
    script:
    """

    vireo \
        -c $sample_outdir \
        -N $params.numsamples \
        -o $sample_outdir/vireo_results \
        --randSeed 2 \
        -p 20


    """
}


// reference based vireo

process vireo_ref {

    label 'process_high'

    publishDir params.publishDir, mode: 'copy'

    container {
    container = 'ghcr.io/cregan727/decoder_ring_vireo:latest'
    }

    
    input:
    path bam_cellsnp
    each path(sample_outdir) 

    output:
    path "*/vireo_results"

    script:
    """
    sample_name=\$(basename $sample_outdir)

    vireo \
        -c "$sample_outdir" \
        -d "$bam_cellsnp" \
        -o "\${sample_name}/vireo_results" \
        --forceLearnGT \
        --randSeed 2 \
        -p 20
    """
}

