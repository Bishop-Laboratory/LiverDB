#!/usr/bin/env nextflow

// CONSTANTS -------------------------------------------------------------------
REF_DIR = "Homo_sapiens.GRCh38.103"
FA_URL = "http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FA_FILENAME = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_URL = "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"
GTF_FILENAME = "Homo_sapiens.GRCh38.103.gtf"

DATA_DIR = "data"

RUN_IDS = file("SRR_list.txt").readLines()
RUN_IDS.removeAll { !it } // remove empty elements from list
run_ch = Channel.value(RUN_IDS)

// PROCESSES -------------------------------------------------------------------
process download_sra {
  
  publishDir "${DATA_DIR}", mode: "move"
  
  input:
  val id from run_ch.flatten()
  
  output:
  file "${id}/${id}.sra" into sra_ch
  
  stub:
  """
  mkdir ${id}
  echo prefetch ${id} > ${id}/${id}.sra
  """
}


process download_assembly {
  
  publishDir "${REF_DIR}", mode: "move"

  output:
  tuple "${FA_FILENAME}", "${GTF_FILENAME}" into assembly_ch
  file "${GTF_FILENAME}" into gtf_align_ch

  stub:
  """
  echo wget -P ${FA_URL} > ${FA_FILENAME}
  echo gunzip ${FA_FILENAME}.gz >> ${FA_FILENAME}
  echo wget -P ${GTF_URL} > ${GTF_FILENAME}
  echo gunzip ${GTF_FILENAME}.gz >> ${GTF_FILENAME}
  """
}


process sra_to_fq {
  
  publishDir "${DATA_DIR}/${sra.simpleName}", mode: "move"
  
  input:
  file sra from sra_ch
  
  output:
  file "${sra}.fastq" into fq_ch
  
  stub:
  """
  ls -lh ${sra} > ${sra}.fastq
  echo fasterq-dump --split-files ${sra} >> ${sra}.fastq
  """
}


process generate_star_index {
  
  publishDir "${REF_DIR}", mode: "move"
  
  input:
  tuple file(fa), file(gtf) from assembly_ch
  
  output:
  path "STAR_index/*" into index_ch
  
  stub:
  """
  mkdir STAR_index
  echo \
  STAR --runMode genomeGenerate \
       --runThreadN 8 \
       --genomeDir STAR_index \
       --genomeFastaFiles ${fa} \
       --sjdbGTFfile ${gtf} > STAR_index/SA
  cp STAR_index/SA STAR_index/SAindex
  cp STAR_index/SA STAR_index/Genome
  """
}


process fastp {
  
  publishDir "${DATA_DIR}/${fq.simpleName}", mode: "move"
  
  input:
  file fq from fq_ch
  
  output:
  file "${fq.simpleName}.trimmed.fastq" into trimmed_fq_ch
  file "${fq.simpleName}-fastp.html"
  file "${fq.simpleName}-fastp.json"
  
  stub:
  """
  ls -lh ${fq} > ${fq.simpleName}.trimmed.fastq
  echo \
  fastp -i ${fq} \
        -o ${fq.simpleName}.trimmed.fastq \
        -h ${fq.simpleName}-fastp.html \
        -j ${fq.simpleName}-fastp.json >> ${fq.simpleName}.trimmed.fastq
  cp ${fq.simpleName}.trimmed.fastq ${fq.simpleName}-fastp.html
  cp ${fq.simpleName}.trimmed.fastq ${fq.simpleName}-fastp.json
  """
}


process star_align_reads {
  
  publishDir "${DATA_DIR}/${tfq.simpleName}", mode: "move"
  
  input:
  path index_files from index_ch
  file gtf from gtf_align_ch
  file tfq from trimmed_fq_ch
  
  output:
  file "${tfq.simpleName}_fakeSTARalignRes"
  
  stub:
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  ls -lh > ${tfq.simpleName}_fakeSTARalignRes
  echo "Now check contents of STAR_index:" >> ${tfq.simpleName}_fakeSTARalignRes
  ls -lh STAR_index >> ${tfq.simpleName}_fakeSTARalignRes
  echo \
  STAR --runMode alignReads --runThreadN 8 --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --readFilesIn ${tfq} \
       --outFileNamePrefix ${tfq.simpleName}_ >> ${tfq.simpleName}_fakeSTARalignRes
  """
}

