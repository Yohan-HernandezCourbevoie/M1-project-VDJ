import os
import glob

configfile: "config.json"

DATA_DIR=config["data_dir"]
ACCESSION_FILES=config['sra_accession']
FLASH_DIR=config["flash_dir"]

SAMPLES=[accession.rstrip() for accession in open(ACCESSION_FILES).readlines()]


rule all:
    input:
        extended=expand(FLASH_DIR+"/{sample}.extendedFrags.fastq.gz", sample=SAMPLES), notcombined=expand(FLASH_DIR+"/{sample}.notCombined_2.fastq.gz", sample=SAMPLES)

def sra_url(srr):
    prefix=srr[:6]
    url=''
    n = len(srr)
    if n == 9:
        url = srr+"/"
    elif n > 9:
        end=srr[9:].zfill(3)
        url=end+"/"+srr+"/"
    return "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"+prefix+"/"+url

rule get_data:
    output:
        sample1=DATA_DIR+"/{accession}_1.fastq.gz", sample2=DATA_DIR+"/{accession}_2.fastq.gz"
    run:
        shell("wget -P {DATA_DIR} -N "+sra_url(wildcards.accession)+"/{wildcards.accession}_{{1,2}}.fastq.gz")

rule run_flash:
    input:
        sample1=DATA_DIR+"/{accession}_1.fastq.gz", sample2=DATA_DIR+"/{accession}_2.fastq.gz"
    output:
        extended=FLASH_DIR+"/{accession}.extendedFrags.fastq.gz", notcombined=FLASH_DIR+"/{accession}.notCombined_2.fastq.gz"
    run:
        shell("flash --suffix=gz -d {FLASH_DIR} {DATA_DIR}/{wildcards.accession}_1.fastq.gz {DATA_DIR}/{wildcards.accession}_2.fastq.gz -o {wildcards.accession}")