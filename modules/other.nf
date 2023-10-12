process process_bam_samplesheet {

    input:
    file params.samplesheet_bams

    output:
    path "bamdir"

    script:
    """

    mkdir bamdir
    # Process the user input samplesheet
    while IFS=',' read -r bam_path bam_name; do
#        ln -s \$bam_path bamdir/\$bam_name.bam
#        ln -s \$bam_path.bai bamdir/\$bam_name.bam.bai
         cp \$bam_path bamdir/\$bam_name.bam 
         cp \$bam_path.bai bamdir/\$bam_name.bam.bai
    done < $params.samplesheet_bams
    # Generate bamsamplesheet.txt

    ls bamdir/*.bam |paste -s -d ',' | sed 's/,\$//' > bamdir/bamsamplesheet.txt
    ls bamdir/*.bam | cut -d "/" -f 2 | sed 's/.bam//' | paste -s -d ',' | sed 's/,\$//' >> bamdir/bamsamplesheet.txt
    """
}


process testing123 {
  cpus 20

  input:
    file 'bamsamplesheet.txt'
    params.regionvcf

  output:
    stdout

  script:
    """
    echo \$(head -n 1 bamsamplesheet.txt)
    
    """
}
