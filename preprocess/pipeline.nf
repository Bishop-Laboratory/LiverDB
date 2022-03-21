#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTANTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

OUTDIR = "raw_counts"
REFS_OUTDIR = "refs"

SKIP_EXISTENT = true  // Skip the processing of a run if its run ID exists in OUTDIR
DELETE_INTERMEDIATES = true  // Delete large intermediate files during workflow

METADATA_CSV = "metadata/metadata.csv"

MAX_RETRIES = 3  // Max number of retries on fail for every process instance

PREFETCH_MAXFORKS   = 2  // Limit parallel prefetch calls so ncbi doesn't get mad
FASTQ_MAXFORKS      = 8  // Limit parallel calls for RAM/IO constraints
FASTP_MAXFORKS      = 8
STAR_ALIGN_MAXFORKS = 6

STAR_GENOME_THREADS = 32
STAR_ALIGN_THREADS  = 16

FA_URL       = "http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FA_FILENAME  = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_URL      = "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"
GTF_FILENAME = "Homo_sapiens.GRCh38.103.gtf"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL INITIALIZATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get run IDs that exist in OUTDIR, if any
def existent_ids = []
fh = new File(OUTDIR)
if(fh.exists()) {
  fh.eachFile {
    existent_ids.add(( it =~ /(SRR\d+)/ )[0][1])
  }
}

// Get sample IDs from metadata.csv
@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv
fh = new File(METADATA_CSV)
def run_ids = parseCsv(fh.getText('utf-8'))

// Populate list of runs to process
def se_runs = []  // IDs of single end runs
def pe_runs = []  // IDs of paired end runs
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

// Starting channel
ch_runs = Channel.value(se_runs + pe_runs)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_FA {
  
  publishDir "${REFS_OUTDIR}", mode: "symlink"
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES

  output:
  file "${FA_FILENAME}" into ch_fa

  script:
  """
  wget ${FA_URL}
  gunzip ${FA_FILENAME}.gz
  """

  stub:
  """
  touch ${FA_FILENAME}
  """
}


process DOWNLOAD_GTF {
  
  publishDir "${REFS_OUTDIR}", mode: "symlink"
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES

  output:
  file "${GTF_FILENAME}" into ch_gtf

  script:
  """
  wget ${GTF_URL}
  gunzip ${GTF_FILENAME}.gz
  """

  stub:
  """
  touch ${GTF_FILENAME}
  """
}


process STAR_GENOME_GENERATE {
  
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file fa from ch_fa
  file gtf from ch_gtf
  
  output:
  path "STAR_index/*" into ch_index
  
  script:
  """
  STAR --runMode genomeGenerate \
       --runThreadN ${STAR_GENOME_THREADS} \
       --genomeDir STAR_index \
       --genomeFastaFiles ${fa} \
       --sjdbGTFfile ${gtf}  
  """
  
  stub:
  """
  mkdir STAR_index
  ls > STAR_index/SA
  cp STAR_index/SA STAR_index/SAindex
  cp STAR_index/SA STAR_index/Genome
  """
}


process PREFETCH {
  
  maxForks PREFETCH_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  val id from ch_runs.flatten()
  
  output:
  file "${id}/${id}.sra" into ch_sra
  
  script:
  """
  prefetch ${id}
  """
  
  stub:
  """
  mkdir ${id}
  touch ${id}/${id}.sra
  """
}


process FASTQ {
  
  maxForks FASTQ_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file sra from ch_sra
  
  output:
  file "${id}_{1,2}.fastq" into ch_fq
  
  script:
  id = sra.simpleName
  """
  fastq-dump --split-files ${sra}
  """
  
  stub:
  id = sra.simpleName
  if(se_runs.contains(id))
    """
    ls > ${id}_1.fastq
    """
  else
    """
    ls > ${id}_1.fastq
    cp ${id}_1.fastq ${id}_2.fastq
    """
}


process FASTP {
  
  maxForks FASTP_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  file fq from ch_fq
  
  output:
  file "${id}_{1,2}.trimmed.fastq" into ch_tfq
  
  script:
  if(!(fq instanceof List)) {  // handle single end
    id = fq.simpleName.split("_")[0]
    """
    fastp -i ${fq} \
          -o ${id}_1.trimmed.fastq
    
    if ${DELETE_INTERMEDIATES}; then
      readlink ${fq} | xargs rm --
      rm ${fq}
    fi
    """
  } else {  // handle paired end
    id = fq[0].simpleName.split("_")[0]
    """
    fastp -i ${fq[0]} \
          -I ${fq[1]} \
          -o ${id}_1.trimmed.fastq \
          -O ${id}_2.trimmed.fastq
    
    if ${DELETE_INTERMEDIATES}; then
      readlink ${fq} | xargs rm --
      rm ${fq}
    fi
    """
  }
  
  stub:
  if(!(fq instanceof List)) {  // handle single end
    id = fq.simpleName.split("_")[0]
    """
    ls > ${id}_1.trimmed.fastq
    
    if ${DELETE_INTERMEDIATES}; then
      readlink ${fq} | xargs rm --
      rm ${fq}
    fi
    """
  } else {  // handle paired end
    id = fq[0].simpleName.split("_")[0]
    """
    ls > ${id}_1.trimmed.fastq
    cp ${id}_1.trimmed.fastq ${id}_2.trimmed.fastq
    
    if ${DELETE_INTERMEDIATES}; then
      readlink ${fq} | xargs rm --
      rm ${fq}
    fi
    """
  }
}


process STAR_ALIGN_READS {
  
  publishDir "${OUTDIR}", mode: "symlink"
  
  maxForks STAR_ALIGN_MAXFORKS
  errorStrategy 'retry'
  maxRetries MAX_RETRIES
  
  input:
  path index_files from ch_index
  file gtf from ch_gtf
  file tfq from ch_tfq
  
  output:
  file "${id}_ReadsPerGene.out.tab" into raw_counts
  
  script:
  if(!(tfq instanceof List)) {
    id = tfq.simpleName.split("_")[0]
  } else {
    id = tfq[0].simpleName.split("_")[0]
  }
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
  if(!(tfq instanceof List)) {
    id = tfq.simpleName.split("_")[0]
  } else {
    id = tfq[0].simpleName.split("_")[0]
  }
  """
  mkdir STAR_index
  mv ${index_files} STAR_index
  ls > ${id}_ReadsPerGene.out.tab

  if ${DELETE_INTERMEDIATES}; then
    readlink ${tfq} | xargs rm --
    rm ${tfq}
  fi
  """
}


