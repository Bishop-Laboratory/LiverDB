#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTANTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

OUTDIR = "out"
SKIP_EXISTENT = true  // Skip the processing of a run if its run ID exists in OUTDIR
DELETE_INTERMEDIATES = true  // Delete large intermediate files during workflow

METADATA_CSV = "metadata/metadata.csv"

MAX_RETRIES = 0  // Allow all process instances to retry on fail

STAR_INDEX_THREADS = 16
STAR_ALIGN_THREADS = 16

PREFETCH_SE_MAXFORKS   = 2  // Limit parallel prefetch calls so ncbi doesn't get mad
PREFETCH_PE_MAXFORKS   = 2
FASTQ_SE_MAXFORKS      = 3  // Limit parallel fastq-dumps for RAM/IO constraints
FASTQ_PE_MAXFORKS      = 6
STAR_ALIGN_SE_MAXFORKS = 2  // Limit parallel STAR alignReads for RAM/IO constraints
STAR_ALIGN_PE_MAXFORKS = 4

FA_URL       = "http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FA_FILENAME  = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_URL      = "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"
GTF_FILENAME = "Homo_sapiens.GRCh38.103.gtf"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL INITIALIZATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get run IDs that exist in OUTDIR
def existent_ids = []
fh = new File(OUTDIR)
fh.eachFile {
  existent_ids.add(( it =~ /(SRR\d+)/ )[0][1])
}

// Get sample IDs from metadata.csv
@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv
fh = new File(METADATA_CSV)
def run_ids = parseCsv(fh.getText('utf-8'))

// Populate list of runs to process
def se_runs = []
def pe_runs = []
for(line in run_ids) {
  id = line.sample_id
  if(SKIP_EXISTENT && existent_ids.contains(id)) {
    continue
  } else if(line.paired_end == "FALSE") {
    se_runs.add(id)
  } else {
    pe_runs.add(id)  
  }
}

// Starting channels
ch_se_runs = Channel.value(se_runs)
ch_pe_runs = Channel.value(pe_runs)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_FA {
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES

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
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES

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
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file fa from ch_fa_index
  file gtf from ch_gtf_index
  
  output:
  path "STAR_index/*" into ch_index
  
  script:
  """
  STAR --runMode genomeGenerate \
       --runThreadN ${STAR_INDEX_THREADS} \
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
       --runThreadN ${STAR_INDEX_THREADS} \
       --genomeDir STAR_index \
       --genomeFastaFiles ${fa} \
       --sjdbGTFfile ${gtf} >> STAR_index/SA
  cp STAR_index/SA STAR_index/SAindex
  cp STAR_index/SA STAR_index/Genome
  """
}


process DOWNLOAD_SRA_SE {
  
  maxForks PREFETCH_SE_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
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
  
  maxForks PREFETCH_PE_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
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


process FASTQ_SE {
  
  maxForks FASTQ_SE_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file sra from ch_sra_se
  
  output:
  file "${sra.simpleName}_1.fastq" into ch_fq_se
  
  script:
  """
  fastq-dump --split-files ${sra}
  """
  
  stub:
  """
  ls -lh > ${sra.simpleName}_1.fastq
  echo fastq-dump --split-files ${sra.simpleName}_1 >> ${sra.simpleName}_1.fastq
  """
}


process FASTQ_PE {
  
  maxForks FASTQ_PE_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file sra from ch_sra_pe
  
  output:
  tuple "${sra.simpleName}_1.fastq", "${sra.simpleName}_2.fastq" into ch_fq_pe
  
  script:
  """
  fastq-dump --split-files ${sra}
  """
  
  stub:
  """
  ls -lh > ${sra.simpleName}_1.fastq
  echo fastq-dump --split-files ${sra.simpleName} >> ${sra.simpleName}_1.fastq
  cp ${sra.simpleName}_1.fastq ${sra.simpleName}_2.fastq
  """
}


process FASTP_SE {
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file fq from ch_fq_se
  
  output:
  file "${fq.simpleName}.trimmed.fastq" into ch_tfq_se
  
  script:
  id = fq.simpleName.split("_")[0]
  """
  fastp -i ${fq} \
        -o ${id}.trimmed.fastq

  if ${DELETE_INTERMEDIATES}; then
    readlink ${fq} | xargs rm --
    rm ${fq}  
  fi
  """
  
  stub:
  id = fq.simpleName.split("_")[0]
  """
  ls -lh ${fq} > ${fq.simpleName}.trimmed.fastq
  echo \
  fastp -i ${fq} \
        -o ${id}.trimmed.fastq >> ${fq.simpleName}.trimmed.fastq

  if ${DELETE_INTERMEDIATES}; then
    readlink ${fq} | xargs rm --
    rm ${fq}  
  fi
  """
}


process FASTP_PE {
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  tuple file(fq1), file(fq2) from ch_fq_pe
  
  output:
  tuple "${id}_1.trimmed.fastq", "${id}_2.trimmed.fastq" into ch_tfq_pe
  
  script:
  id = fq1.simpleName.split("_")[0]
  """
  fastp -i ${fq1} \
        -I ${fq2} \
        -o ${id}_1.trimmed.fastq \
        -O ${id}_2.trimmed.fastq
        
  if ${DELETE_INTERMEDIATES}; then
    readlink ${fq1} | xargs rm --
    rm ${fq1}
    readlink ${fq2} | xargs rm --
    rm ${fq2} 
  fi
  """
  
  stub:
  id = fq1.simpleName.split("_")[0]
  """
  ls -lh > ${id}_1.trimmed.fastq
  echo \
  fastp -i ${fq1} \
        -I ${fq2} \
        -o ${id}_1.trimmed.fastq \
        -O ${id}_2.trimmed.fastq  >> ${id}_1.trimmed.fastq
  cp ${id}_1.trimmed.fastq ${id}_2.trimmed.fastq
  
  if ${DELETE_INTERMEDIATES}; then
    readlink ${fq1} | xargs rm --
    rm ${fq1}
    readlink ${fq2} | xargs rm --
    rm ${fq2} 
  fi
  """
}


process STAR_ALIGN_SE {
  
  publishDir "${OUTDIR}", mode: "symlink"
  
  maxForks STAR_ALIGN_SE_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  path index_files from ch_index
  file gtf from ch_gtf_align
  file tfq from ch_tfq_se
  
  output:
  file "${id}_ReadsPerGene.out.tab"
  
  script:
  id = tfq.simpleName.split("_")[0]
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  STAR --runMode alignReads \
       --runThreadN ${STAR_ALIGN_THREADS} \
       --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --readFilesIn ${tfq} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${id}_
  
  if ${DELETE_INTERMEDIATES}; then
    readlink ${tfq} | xargs rm --
    rm ${tfq}
    rm *_Aligned.sortedByCoord.out.bam
  fi
  """
  
  stub:
  id = tfq.simpleName.split("_")[0]
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  ls -lh > ${id}_ReadsPerGene.out.tab
  echo "Now check contents of STAR_index:" >> ${id}_ReadsPerGene.out.tab
  ls -lh STAR_index >> ${id}_ReadsPerGene.out.tab
  echo \
  STAR --runMode alignReads \
       --runThreadN ${STAR_ALIGN_THREADS} \
       --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --readFilesIn ${tfq} \
       --outFileNamePrefix ${id}_ >> ${id}_ReadsPerGene.out.tab
       
  if ${DELETE_INTERMEDIATES}; then
    readlink ${tfq} | xargs rm --
    rm ${tfq}
  fi
  """
}


process STAR_ALIGN_PE {
  
  publishDir "${OUTDIR}", mode: "symlink"
  
  maxForks STAR_ALIGN_PE_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  path index_files from ch_index
  file gtf from ch_gtf_align
  tuple file(tfq1), file(tfq2) from ch_tfq_pe
  
  output:
  file "${id}_ReadsPerGene.out.tab"
  
  script:
  id = tfq1.simpleName.split('_')[0]
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  STAR --runMode alignReads \
       --runThreadN ${STAR_ALIGN_THREADS} \
       --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --readFilesIn ${tfq1} ${tfq2} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${id}_
       
  if ${DELETE_INTERMEDIATES}; then
    readlink ${tfq1} | xargs rm --
    readlink ${tfq2} | xargs rm --
    rm ${tfq1} ${tfq2}
    rm *_Aligned.sortedByCoord.out.bam
  fi
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
  STAR --runMode alignReads \
       --runThreadN ${STAR_ALIGN_THREADS} \
       --quantMode GeneCounts \
       --genomeDir STAR_index \
       --sjdbGTFfile ${gtf} \
       --readFilesIn ${tfq1} ${tfq2} \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${id}_ >> ${id}_ReadsPerGene.out.tab
       
  if ${DELETE_INTERMEDIATES}; then
    readlink ${tfq1} | xargs rm --
    readlink ${tfq2} | xargs rm --
    rm ${tfq1} ${tfq2}
  fi
  """
}




