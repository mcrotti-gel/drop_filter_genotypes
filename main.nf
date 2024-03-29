#!/usr/bin/env nextflow

/* 
 * Set up variables
 */

Channel
      .fromPath(params.vcf_list)
      .ifEmpty { exit 1, "Cannot find input file : ${params.vcf_list}" }
      .splitCsv(skip:1)
      .map { row -> tuple(row[0], file(row[1]), file(row[2]))}
      .take( params.number_of_files_to_process )
      .set { vcf_input }

/*
 * Start pipeline
 */

process drop_genotypes {

  tag "$sampleID"

   input:
   tuple val(sampleID), file(vcf), file(index) from vcf_input

   output:
   tuple val(sampleID), file("*_sites.vcf.gz"), file("*_sites.vcf.gz.csi") into drop_geno_ch

   script:

   """
   bcftools view -G -O z -o ${sampleID}_sites.vcf.gz ${vcf} 
   bcftools index ${sampleID}_sites.vcf.gz
   """

}

process medianGQ_filter {

  tag "sampleID"
  publishDir "${params.outdir}/subset", mode: 'copy'
  
  input:
  tuple val(sampleID), file(vcf), file(index) from drop_geno_ch

  output:
  file("*_medianGQ.tsv") into medianGQ_filtered_ch

  script:

  """
  bcftools query -i 'medianGQ > 30' -f '%CHROM\t%POS\t%REF\t%ALT\t%medianGQ\n' ${vcf} > ${sampleID}_medianGQ.tsv
  """
}
