# How does nf-decoder-ring work?

nf-decoder-ring is a [Nextflow](https://www.nextflow.io/) based pipeline
that takes advantage of docker containers to make running a [CellSNP](https://cellsnp-lite.readthedocs.io/en/latest/manual.html) -\>
[Vireo](https://vireosnp.readthedocs.io/en/latest/manual.html) genetic demultiplexing workflow simple and easy. The inputs are a
sample sheet containing the paths to the 10x bam files we want to
demultiplex and if available a second sample sheet defining the paths to
any reference bam files and metadata. The pipeline is designed to work
for 3 common use cases: reference free demultiplexing, reference based
demultiplexing when the reference includes 1 bam per sample, reference
based demultiplexing when the reference is from a bam where each sample
is barcoded.

![nf-decoder-ring-wf](https://github.com/cregan727/nf-decoder-ring/assets/68451521/7486fd6b-682a-46d7-aa70-9320cde25daf)

Below we'll walk through an example of an experiment that
uses plate based RNA-seq as a reference and how nf-decoder-ring can be
used to demultiplex the samples.

## Plate Based RNA-seq example

Alongside a mixed donor 10x Genomics 5' single cell experiment, we
pelleted 1 Million cells from each donor, trizol extracted RNA and
performed barcoded low input RNA-seq to use as a reference. We used
barcoded RT primers to keep the low input RNA-seq costs down as it
allows the samples to be pooled immediately after the RT instead of
processing them as independent samples. After the library was prepared
it was sequenced with the scRNA-seq libraries we want to decode on a
NextSeq 2000 P3 to \~125 Million reads total.

Barcoded RT primer design Design:

```         
    Illumina R1          BC          UMI
CTACACGACGCTCTTCCGATCT-AACGTGAT-NNNNNNNNNNVVVVVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN
```

In our experiment these were the barcodes assigned to our samples:

AACAACCA,Celiac

AACCGAGA,Healthy

AACGCTTA,Osteo

AAGACGGA,Rheuma


## Pre-processing ahead of nf-decoder-ring

nf-decoder-ring expects bam files as input for both the single cell
libraries and the reference. In this example we first had to map the low
input RNA-seq data to the human genome using STARsolo which will
generate the barcoded bam file we need to be the reference for
nf-decoder-ring.

GENEREAD = /path/to/R2

CBCUMIREAD = /path/to/R1

WHITELIST = a list of the potential barcodes in the low input RNA-seq

BARCODELENGTH = 8

UMISTART = 9

UMILEN = 16

R1LENGTH = in our case 26 since these were sequenced with the 10x 5'
libraries that required a 26bp R1

```         
STAR \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN=20 \
        --genomeDir $GENOME \
        --readFilesCommand zcat \
        --readFilesIn $GENEREAD $CBCUMIREAD \
        --soloType Droplet \
        --soloCBwhitelist $WHITELIST \
        --soloCBstart 1 \
        --soloCBlen $BARCODELENGTH \
        --soloUMIstart $UMISTART \
        --soloUMIlen 16 \ #$UMILEN \
        --soloBarcodeReadLength $R1LENGTH \
```

## What happens when you run nf-decoder-ring?

nf-decoder-ring runs as a single command, the user has to provide the
following:

\--samplesheet - path to a sample sheet containing samples names and the
paths to their 10x bams

\--regionvcf - path to a vcf of candidate SNPs to use for the CellSNP
genotyping

\--publishDir - path to the directory you want the pipeline outputs to
end up in

\--workflow - which workflow you want to use (options are "ref_free",
"ref_bams", or "ref_plate")

\--samplesheet_plate - in order to use nf-decoder-ring provide a
samplesheet with path to the reference bam and metadata

-profile - use either singularity or docker to run the pipeline without
a local installation of cellsnp and vireo

-resume - use this option to allow you to resume a previous running of
the pipeline (if something fails you have the opportunity to fix it and
rerun it from where it stopped instead of rerunning the whole thing)

```         
nextflow run cregan727/nf-decoder-ring -r main --samplesheet $PWD/samples.analysis.csv \
        --regionvcf "path/to/vcf.gz"   \
        --publishDir $PWD/published_results \
        --workflow "ref_plate" \
        --samplesheet_plate $PWD/jplate_samplesheet.csv \
        -profile singularity -resume 
```

The first thing that happens after you run your command is Nextflow will
parse your sample sheets and begin the ref_plate workflow. First it
splits the reference sample sheet into columns which correspond to the
path to the bam file and the path to a sample metadata file. It uses this and the regionvcf
file as input into the cellsnp_lite_plate process.

It also similarly parses the samplesheet from the 10x samples to
run the process cellsnp_lite_10x.

When both of these processes are complete it runs vireo_ref using the
outputs from both processes.

```         
input_data = Channel.fromPath(params.samplesheet)
                      .splitCsv(header: true, sep: ',')
                      
input_ref = Channel.fromPath(params.samplesheet_plate)
                      .splitCsv(header: false, sep: ',')

workflow ref_plate {


    // take input as bam files and run cellsnp-lite
                            
    bam_cellsnp = input_ref.map { sheet ->
    tuple(file(sheet[0]), file(sheet[0] + ".bai"), file(sheet[1]), file(params.regionvcf))
} | cellsnp_lite_plate


    //perform cellsnp-lite on the 10x bams provided in the samplesheet
    
    wf1_out = input_data.map { row ->
    tuple(row['Sample'], file(row['bam']), file(row['bam'] + ".bai"), file(row['barcode']), file(params.regionvcf))
} | cellsnp_lite_10x


    //run reference based vireo
    vireo_ref(bam_cellsnp, wf1_out)

}
```

## CellSNP

Here is an example of what the cellsnp_lite_10x process looks like.

At the top it defines the publish directory, the number of CPUs and
where I want it to pull its' docker container from. Then I specify the
inputs which are the sample name, the bam, bam index, cellranger called
cell barcodes, and the region vcf file. Then I define the outputs I want
to capture as the path "\${sample}\_outdir/". Finally I define the
actual script to run cellsnp-lite which is a command line tool.

cellSNP genotypes either the single cells or the reference and generates
a vcf file/mtx files for use in vireo.

```         
process cellsnp_lite_10x {

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
```

## Vireo

Vireo takes the cellsnp called genotypes and the references and assigns
each cell to a sample. Here is what the vireo command looks like:

The inputs are the reference vcf (here named bam_cellsnp since that's the process it comes from)
and each path(sample_outdir). The 'each' indicates that I want to run this
process for each item in that channel which in this case is a 1 cellsnp output directory for each 10x lane. 
The output is the path to the vireo output directory.

```         
process vireo_ref {

    publishDir params.publishDir, mode: 'copy'

    container {
    container = 'swarbricklab/vireo_snp:0.5.6'
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
        --randSeed 2 \
        -p 20
    """
}
```

## Results

The results from both CellSNP and Vireo are published to the user
defined Publish Directory. The CellSNP directory contains the
intermediate files used by Vireo and the Vireo output directory contains
files about your sample assignment.

Summary File:

```         
Var1    Freq
Celiac  1745
Healthy 2406
Osteo   2425
Rheuma  1972
doublet 647
unassigned      994
```

A PDF heatmap of how genetically similar the donors are to each other:
<img width="819" alt="GenoProbDelta" src="https://github.com/cregan727/nf-decoder-ring/assets/68451521/286434c0-8c61-4ba5-b60b-029bd6e84bf3">


Two TSVs of the probability each cell is each individual donor or a doublet of two donors

A log file

And the donor assignments for each cell in a file like this:

```         
cell    donor_id        prob_max        prob_doublet    n_vars  best_singlet    best_doublet    doublet_logLikRatio
AAACCTGAGACCTAGG-1      Osteo   1.00e+00        4.90e-12        143     Osteo   Celiac,Osteo    -23.459
AAACCTGAGACTGTAA-1      Celiac  1.00e+00        4.22e-05        58      Celiac  Celiac,Osteo    -7.491
AAACCTGAGGTCATCT-1      Osteo   9.79e-01        1.81e-02        33      Osteo   Celiac,Osteo    -1.408
AAACCTGAGTATTGGA-1      Healthy 1.00e+00        1.50e-05        63      Healthy Healthy,Osteo   -8.527
```

The CellSNP output is also visible in this folder though it is unlikely you'll need anything from this for downstream analysis. Here is an excerpt from the CellSNP documentation explaining what the output files are. 

```
cellSNP.base.vcf.gz: a VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT).

cellSNP.samples.tsv: a TSV file listing cell barcodes or sample IDs.

cellSNP.tag.AD.mtx: a file in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles.

cellSNP.tag.DP.mtx: a file in “Matrix Market exchange formats”, containing the sum of allele depths of the reference and alternative alleles (REF + ALT).

cellSNP.tag.OTH.mtx: a file in “Matrix Market exchange formats”, containing the sum of allele depths of all the alleles other than REF and ALT.
```

When running CellSNP on the reference bam we also generate a VCF file which contains the genotype for each of the samples which is used as an input for Vireo. This can be found in the output folder of the publishDIR.
