#!/usr/bin/env python3
import pandas
from snakemake.remote.HTTP import RemoteProvider as HTTPprovider
from tempfile import mkdtemp

#############
# CONTAINERS #
#############

bbmap = "docker://quay.io/biocontainers/bbmap:38.98--h5c4e2a8_1"
star = "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"

#############
# FUNCTIONS #
#############

#TODO ask what this does?
def maketempdir():
    return Path(mkdtemp(), 'tmp').resolve().as_posix()

###########
# GLOBALS #
###########
sample_table_loc = "data/walsh_sample_table.txt"
reads_dir = 'data/walsh_raw_mRNA_reads'

#genomes downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/
ref_gff_url = ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz')
ref_gff = 'GCF_000001405.40_GRCh38.p14_genomic.gff'

ref_fna_url = ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz')  
ref_fna = 'GCF_000001405.40_GRCh38.p14_genomic.fna'  



#########
# MAIN #
#########
HTTP = HTTPprovider()

sample_table = pandas.read_csv(
    sample_table_loc,
    index_col="sample ID", 
    sep='\t', 
    lineterminator='\r')

#########
# RULES #
#########

rule target:
    input:
        expand('output/star/pass2/{sample}.ReadsPerGene.out.tab',
               sample=paired_sample_names)
        
        
rule star_second_pass:
    input:
        r1 = 'output/trim/{sample}_1.fastq.gz',
        r2 = 'output/trim/{sample}_2.fastq.gz',
        star_reference = 'output/star/star-index',
        junctions = expand('output/star/pass1/{sample}.SJ.out.tab',
                           sample=paired_sample_names)
    output:
        counts = 'output/star/pass2/{sample}.ReadsPerGene.out.tab'
    threads:
        10
    params:
        prefix = 'output/star/pass2/{sample}.'
    log:
        'output/logs/star_second_pass.{sample}.log'
    resources:
        time = 99,
        mem_mb = 64 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {input.star_reference} ' 
        '--sjdbFileChrStartEnd {input.junctions} ' 
        '--outSAMtype None ' 
        '--quantMode GeneCounts ' 
        '--readFilesCommand zcat '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '--outTmpDir ' + maketempdir() + ' '
        '&> {log}'
        

        
      
rule star_first_pass:
    input:
        r1 = 'output/trim/{sample}_1.fastq.gz',
        r2 = 'output/trim/{sample}_2.fastq.gz',
        star_reference = 'output/star/star-index' 
    output:
        sjdb = 'output/star/pass1/{sample}.SJ.out.tab'
    threads:
        10
    params:
        prefix = 'output/star/pass1/{sample}.'  
    log:
        'output/logs/star_first_pass.{sample}.log'
    resources:
        time = 99,
        mem_mb = 64 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {input.star_reference} '
        '--outSJfilterReads Unique ' 
        '--outSAMtype None '
        '--readFilesCommand zcat '
        '--readFilesIn {input.r1} {input.r2} ' 
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'
        
        
rule star_index:
    input:
        gff = f'output/ref/{ref_gff}', 
        fasta = f'output/ref/{ref_fna}'
    output:
        directory('output/star/star-index')
    params:
        outdir = 'output/star/star-index'
    log:
        'output/logs/star_index.log'
    threads:
        10
    resources:
        time = 99,
        mem_mb = 64 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gff} ' 
        '--genomeSAindexNbases 12 '
        '--outTmpDir ' + maketempdir() + ' '
        '--sjdbGTFtagExonParentTranscript Parent '
        '--sjdbGTFtagExonParentGene gene '
        '&> {log}'
        
rule download_ref_fna:
    input:
        HTTP.remote(ref_fna_url, keep_local=True)
    output:
        f'output/ref/{ref_fna}'
    log:
        'output/logs/ref_fna.log'
    shell:
        'gunzip -c {input} > {output} '
    
rule download_ref_gff:
    input:
        HTTP.remote(ref_gff_url, keep_local=True)
    output:
        f'output/ref/{ref_gff}'
    log:
        'output/logs/ref_gff.log'
    shell:
        'gunzip -c {input} > {output} '

rule trim:
    input:
        r1 = Path(reads_dir, "{sample}_1.fastq.gz").resolve(),
        r2 = Path(reads_dir, "{sample}_2.fastq.gz").resolve()
        
    output:
        r1 = 'output/trim/{sample}_1.fastq.gz',
        r2 = 'output/trim/{sample}_2.fastq.gz'
        
    log:
        'output/logs/trim.{sample}.log'
    threads:
        10
    resources:
        time = 59,
        mem_mb = 10 * 1000
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'zl=9 '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'
