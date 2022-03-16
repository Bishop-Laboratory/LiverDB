#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTANTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

OUTDIR = "out"

FA_URL       = "http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FA_FILENAME  = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_URL      = "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"
GTF_FILENAME = "Homo_sapiens.GRCh38.103.gtf"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL INITIALIZATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get sample IDs from metadata.csv
@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv
 
fh = new File("metadata/metadata.csv")
def data = parseCsv(fh.getText('utf-8'))

def se_runs = []
def pe_runs = []

for(line in data) {
  if(line.paired_end == "FALSE") {
    se_runs.add(line.sample_id)
  } else {
    pe_runs.add(line.sample_id)  
  }
}

// Starting channels
ch_se_runs   = Channel.value(se_runs)
ch_pe_runs   = Channel.value(pe_runs)
ch_rscript   = Channel.fromPath("scripts/downstream.R")
ch_metadata  = Channel.fromPath("metadata/metadata.csv")
ch_contrasts = Channel.fromPath("metadata/contrasts.csv")



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_FA {

  output:
  file "${FA_FILENAME}" into ch_fa_index

  script:
  """
  wget ${FA_URL}
  gunzip ${FA_FILENAME}.gz
  """

  stub:
  """
  echo wget ${FA_URL} > ${FA_FILENAME}
  echo gunzip ${FA_FILENAME}.gz >> ${FA_FILENAME}
  """
}


process DOWNLOAD_GTF {

  output:
  file "${GTF_FILENAME}" into ch_gtf_index, ch_gtf_align

  script:
  """
  wget ${GTF_URL}
  gunzip ${GTF_FILENAME}.gz
  """

  stub:
  """
  echo wget ${GTF_URL} > ${GTF_FILENAME}
  echo gunzip ${GTF_FILENAME}.gz >> ${GTF_FILENAME}
  """
}


process GENERATE_STAR_INDEX {
  
  input:
  file fa from ch_fa_index
  file gtf from ch_gtf_index
  
  output:
  path "STAR_index/*" into ch_index
  
  script:
  """
  STAR --runMode genomeGenerate \
       --runThreadN 8 \
       --genomeDir STAR_index \
       --genomeFastaFiles ${fa} \
       --sjdbGTFfile ${gtf}  
  """
  
  stub:
  """
  mkdir STAR_index
  ls -lh > STAR_index/SA 
  echo \
  STAR --runMode genomeGenerate \
       --runThreadN 8 \
       --genomeDir STAR_index \
       --genomeFastaFiles ${fa} \
       --sjdbGTFfile ${gtf} >> STAR_index/SA
  cp STAR_index/SA STAR_index/SAindex
  cp STAR_index/SA STAR_index/Genome
  """
}


process DOWNLOAD_SRA_SE {
  
  input:
  val id from ch_se_runs.flatten()
  
  output:
  file "${id}/${id}.sra" into ch_sra_se
  
  script:
  """
  prefetch ${id}
  """
  
  stub:
  """
  mkdir ${id}
  echo prefetch ${id} > ${id}/${id}.sra
  """
}


process DOWNLOAD_SRA_PE {
  
  input:
  val id from ch_pe_runs.flatten()
  
  output:
  file "${id}/${id}.sra" into ch_sra_pe
  
  script:
  """
  prefetch ${id}
  """
  
  stub:
  """
  mkdir ${id}
  echo prefetch ${id} > ${id}/${id}.sra
  """
}


process SRA_TO_FQ_SE {
  
  input:
  file sra from ch_sra_se
  
  output:
  file "${sra}.fastq" into ch_fq_se
  
  script:
  """
  fasterq-dump --split-files ${sra}
  """
  
  stub:
  """
  ls -lh > ${sra}.fastq
  echo fasterq-dump --split-files ${sra} >> ${sra}.fastq
  """
}


process SRA_TO_FQ_PE {
  
  input:
  file sra from ch_sra_pe
  
  output:
  tuple "${sra}_1.fastq", "${sra}_2.fastq" into ch_fq_pe
  
  script:
  """
  fasterq-dump --split-files ${sra}
  """
  
  stub:
  """
  ls -lh > ${sra}_1.fastq
  echo fasterq-dump --split-files ${sra} >> ${sra}_1.fastq
  cp ${sra}_1.fastq ${sra}_2.fastq
  """
}


process FASTP_SE {
  
  input:
  file fq from ch_fq_se
  
  output:
  file "${fq.simpleName}.trimmed.fastq" into ch_tfq_se
  file "${fq.simpleName}-fastp.html"
  file "${fq.simpleName}-fastp.json"
  
  script:
  """
  fastp -i ${fq} \
        -o ${fq.simpleName}.trimmed.fastq \
        -h ${fq.simpleName}-fastp.html \
        -j ${fq.simpleName}-fastp.json
  """
  
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


process FASTP_PE {
  
  input:
  tuple file(fq1), file(fq2) from ch_fq_pe
  
  output:
  tuple "${id}_1.trimmed.fastq", "${id}_2.trimmed.fastq" into ch_tfq_pe
  file "${id}-fastp.html"
  file "${id}-fastp.json"
  
  script:
  id = fq1.simpleName
  """
  fastp -i ${fq1} \
        -I ${fq2} \
        -o ${id}_1.trimmed.fastq \
        -O ${id}_2.trimmed.fastq \
        -h ${id}-fastp.html \
        -j ${id}}-fastp.json
  """
  
  stub:
  id = fq1.simpleName
  """
  ls -lh > ${id}_1.trimmed.fastq
  echo \
  fastp -i ${fq1} \
        -I ${fq2} \
        -o ${id}_1.trimmed.fastq \
        -O ${id}_2.trimmed.fastq \
        -h ${id}-fastp.html \
        -j ${id}}-fastp.json >> ${id}_1.trimmed.fastq
  cp ${id}_1.trimmed.fastq ${id}_2.trimmed.fastq
  cp ${id}_1.trimmed.fastq ${id}-fastp.html
  cp ${id}_1.trimmed.fastq ${id}-fastp.json
  """
}


process STAR_ALIGN_SE {
  
  publishDir "${OUTDIR}", mode: "symlink"
  
  input:
  path index_files from ch_index
  file gtf from ch_gtf_align
  file tfq from ch_tfq_se
  
  output:
  file "${id}_ReadsPerGene.out.tab" into ch_raw_reads_se
  
  script:
  id = tfq.simpleName
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  STAR --runMode alignReads --runThreadN 8 --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --readFilesIn ${tfq} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${id}_
  """
  
  stub:
  id = tfq.simpleName
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  ls -lh > ${id}_ReadsPerGene.out.tab
  echo "Now check contents of STAR_index:" >> ${id}_ReadsPerGene.out.tab
  ls -lh STAR_index >> ${id}_ReadsPerGene.out.tab
  echo \
  STAR --runMode alignReads --runThreadN 8 --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --readFilesIn ${tfq} \
       --outFileNamePrefix ${id}_ >> ${id}_ReadsPerGene.out.tab
  """
}


process STAR_ALIGN_PE {
  
  publishDir "${OUTDIR}", mode: "symlink"
  
  input:
  path index_files from ch_index
  file gtf from ch_gtf_align
  tuple file(tfq1), file(tfq2) from ch_tfq_pe
  
  output:
  file "${id}_ReadsPerGene.out.tab" into ch_raw_reads_pe
  
  script:
  id = tfq1.simpleName.split('_')[0]
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  STAR --runMode alignReads --runThreadN 8 --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --readFilesIn ${tfq1} ${tfq2} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${id}_
  """
  
  stub:
  id = tfq1.simpleName.split('_')[0]
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  ls -lh > ${id}_ReadsPerGene.out.tab
  echo "Now check contents of STAR_index:" >> ${id}_ReadsPerGene.out.tab
  ls -lh STAR_index >> ${id}_ReadsPerGene.out.tab
  echo \
  STAR --runMode alignReads --runThreadN 8 --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --readFilesIn ${tfq1} ${tfq2} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${id}_ >> ${id}_ReadsPerGene.out.tab
  """
}




