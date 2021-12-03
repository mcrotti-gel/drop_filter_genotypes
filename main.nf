#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* 
 * Set up variables
 */


Channel
      .fromPath(params.vcf_list)
      .ifEmpty { exit 1, "Cannot find input file : ${params.vcf_list}" }
      .splitCsv()
      .map { row -> tuple(row[0], file(row[1]), file(row[2]))}
      .take( params.number_of_files_to_process )
      .set { vcf_input }

workflow {
  drop_genotypes( vcf_input )
  AF_filter( drop_genotypes.out )

}

process drop_genotypes {
  
  tag "$sampleID"

  input:
  tuple val(sampleID), file(vcf), file(index)

  output:
  tuple val(sampleID), file("*_sites.vcf.gz"), file("*_sites.vcf.gz.csi")

  script:

   """
   bcftools view -G -O z -o ${sampleID}_sites.vcf.gz ${vcf} 
   bcftools index ${sampleID}_sites.vcf.gz
   """

}

process minGQ_filter {

  tag "sampleID"
  
  input:
  tuple val(sampleID), file(vcf), file(index)

  output:
  file("*_maf02.tsv") 

  script:

  """
  bcftools query -i 'minGQ < 30' -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' ${vcf} > ${sampleID}_minGQ.tsv
  """
}

