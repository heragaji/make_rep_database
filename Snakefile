import sys
import subprocess
import os

configfile:
    "config.json"

subprocess.call("seqkit sort -lr {reads_data} | seqkit head -n {top} > {output}/top.fa".format(reads_data=config['data'],top=config['top'],output=config['output']),shell = True)
top_id = subprocess.check_output('seqkit seq -ni {output}/top.fa'.format(output=config['output']),shell = True).decode(encoding='utf-8').split()
for id in top_id:
    os.makedirs(config['output']+"/"+id, exist_ok=True)

rule all:
    input:
        config['output']+"/cluster_size.csv"

rule minialign:
    input:
        top = config['output']+"/top.fa"
    output:
        config['output']+"/top_vs_reads_sorted.bam"
    params:
        data = config['data']
    shell:
        """
        /home/ozawat/local/src/minialign/minialign -x pacbio -f 0 -m 0.00001 {input.top} {params.data}  -t16 | samtools view -Sb -F 4 | samtools sort > {output}
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
        mkfifo {params.dir}/{wildcards.ref}/hoge
        samtools view -h {input.bam} > {params.dir}/{wildcards.ref}/hoge &
        STACK_YAML={params.realigner}/stack.yaml stack exec -- realigner -i {params.it} {wildcards.ref} < {params.dir}/{wildcards.ref}/hoge | samtools view -Sb | samtools sort >  {params.dir}/{wildcards.ref}/{wildcards.ref}_vs_reads_realigned.bam
        rm {params.dir}/{wildcards.ref}/hoge
        """


rule dump_consensus:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_realigned.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_realigned_count.tsv"
    params:
        dump_dir = config['dump_consensus']
    shell:
        """
        samtools view -h {input} | STACK_YAML={params.dump_dir}/stack.yaml stack exec -- dump-consensus {wildcards.ref} -c > {output}
        """

rule make_consensus:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_realigned_count.tsv"
    output:
        config['output']+"/{ref}/consensus_{ref}_vs_reads_realigned.fa"
    params:
        dir = config['output'],
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/count_to_consensus.py -s {wildcards.ref} -f {input} -g True > {output}
        """

rule coverage:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_realigned.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_depth.txt"
    params:
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/bam-alignment_coverage.py {input} > {output}
        """

rule region:
     input:
         config['output']+"/{ref}/{ref}_vs_reads_depth.txt"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_region.txt"
     params:
         cov = config['coverage'],
         src_dir = config['src_dir']
     shell:
         """
         python {params.src_dir}/depth-high_coverage.py {input}  {params.cov} > {output}
         """

rule terminal:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_realigned.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_terminal.tsv"
    params:
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/bam-alignment_terminal.py {input} > {output}
        """

rule peak:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_terminal.tsv"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_peak.txt"
    params:
        interval = config['interval'],
        peak = config['peak'],
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/terminal-find_peaks.py {input} {params.interval} {params.peak} > {output}
        """

rule cut:
    input:
        region = config['output']+"/{ref}/{ref}_vs_reads_region.txt",
        peak = config['output']+"/{ref}/{ref}_vs_reads_peak.txt"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_cut.bed"
    params:
        cut = config['cut'],
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/peak-cut_region.py {input.region} {input.peak} {params.cut} > {output}
        """

rule gapped_result:
     input:
         bed = config['output']+"/{ref}/{ref}_vs_reads_cut.bed",
         consensus = config['output']+"/{ref}/consensus_{ref}_vs_reads_realigned.fa"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_result_gapped.fa"
     shell:
         """
         seqkit subseq --bed {input.bed} {input.consensus} > {output}
         """

rule each_result:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_result_gapped.fa"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_result.fa"
    params:
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/remove_gap_from_fasta.py {input} > {output}
        """

rule repeat:
     input:
         expand(config['output']+"/{ref}/{ref}_vs_reads_result_gapped.fa",ref=top_id)
     output:
         config['output']+"/repeat.fa"
     shell:
         """
         touch {output}
         cat {input} >> {output}
         """

rule cluster:
    input:
        config['output']+"/repeat.fa"
    output:
        config['output']+"/cluster_size.csv"
    params:
        src_dir = config['src_dir']
    shell:
        """
        python {params.src_dir}/cal_cluster_size.py {input} > {output}
        """
