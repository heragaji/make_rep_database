import sys
import subprocess
import os

configfile:
    "config.json"


subprocess.call("seqkit sort -lr {reads_data} | seqkit head -n {top} > {output}/top.fa".format(reads_data=config['data'],top=config['top'],output=config['output']),shell = True)
top_id = subprocess.check_output('seqkit seq -ni {output}/top.fa'.format(output=config['output']),shell = True).decode(encoding='utf-8').split()
for id in top_id:
    os.makedirs(config['output']+"/"+id, exist_ok=True)
src_dir = srcdir("src")
include_dir = srcdir("after_realigner")

include: include_dir + "/Snakefile"

rule all:
    input:
        config['output']+"/read_masked.fa",
        config['output']+"/weak.paf",
        config['output']+"/rep_vs_rep.paf"


rule minialign:
    input:
        top = config['output']+"/top.fa"
    output:
        config['output']+"/top_vs_reads_sorted.bam"
    params:
        data = config['data']
    threads: 16
    shell:
        """
        minialign -x pacbio -f 0 -m 0.00001 {input.top} {params.data}  -t{threads} | samtools view -Sb -F 4 | samtools sort > {output}
        """

rule realigner:
    input:
        bam = config['output']+"/top_vs_reads_sorted.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_realigned.bam"
    params:
        realigner = config['realigner'],
        it = config['iteration'],
        dir = config['output'],
        top_fa = config['output']+"/top.fa"

    shell:
        """
        if test -e {params.dir}/{wildcards.ref}/hoge; then
            rm {params.dir}/{wildcards.ref}/hoge
        fi
        mkfifo {params.dir}/{wildcards.ref}/hoge
        samtools view -h {input.bam} > {params.dir}/{wildcards.ref}/hoge &
        STACK_YAML={params.realigner}/stack.yaml stack exec -- realigner -i {params.it} {wildcards.ref} < {params.dir}/{wildcards.ref}/hoge | samtools view -Sb | samtools sort >  {params.dir}/{wildcards.ref}/{wildcards.ref}_vs_reads_realigned.bam
        rm {params.dir}/{wildcards.ref}/hoge
        """
