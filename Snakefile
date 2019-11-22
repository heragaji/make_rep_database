import sys
import subprocess
import os

configfile:
    "config.json"

subprocess.call("seqkit sort -lr {reads_data} | seqkit head -n {top} > {output}/top{top}.fa".format(reads_data=config['data'],top=config['top'],output=config['output']),shell = True)
top_id = subprocess.check_output('seqkit seq -ni {output}/top{top}.fa'.format(output=config['output'],top=config['top']),shell = True).decode(encoding='utf-8').split()
for id in top_id:
    os.makedirs(config['output']+"/"+id, exist_ok=True)
src_dir = srcdir("src")

rule all:
    input:
        config['output']+"/repeat_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"

# rule minialign:
#     output:
#         config['output']+"/top_vs_reads_sorted_top"+config['top']+".bam"
#     params:
#         data = config['data']
#     threads: 16
#     shell:
#         """
#         minialign -x pacbio -f 0 -m 0.00001 {config[output]}/top.fa {params.data}  -t{threads} | samtools view -Sb -F 4 | samtools sort > {output}
#         """
rule split_fasta:
    input:
        config['output']+"/top"+config['top']+".fa"
    output:
        config['output']+"/split"
    shell:
        """
        python {src_dir}/split_fasta.py {input} {config[output]} {config[top]} > {output}
        """

rule minialign:
    input:
        tmp=config['output']+"/split",
    output:
        config['output']+"/{ref}/top_vs_reads_sorted_top"+config['top']+"_{ref}.bam"
    params:
        data = config['data'],
        ref = config['output']+"/{ref}/top"+config['top']+"_{ref}.fa"
    threads: 16
    shell:
        """
        minialign -x pacbio -f 0 -m 0.00001 {params.ref} {params.data}  -t{threads} | samtools view -Sb -F 4 | samtools sort > {output}
        """

# rule realigner:
#     input:
#         config['output']+"/{ref}/top_vs_reads_sorted_top"+config['top']+"_{ref}.bam"
#     output:
#         config['output']+"/{ref}/{ref}_vs_reads_realigned_top"+config['top']+"_it"+config['iteration']+".bam"
#     params:
#         realigner = config['realigner'],
#         it = config['iteration'],
#         dir = config['output'],
#     shell:
#         """
#         if test -e {params.dir}/{wildcards.ref}/hoge; then
#             rm {params.dir}/{wildcards.ref}/hoge
#         fi
#         mkfifo {params.dir}/{wildcards.ref}/hoge
#         samtools view -h {input} > {params.dir}/{wildcards.ref}/hoge &
#         STACK_YAML={params.realigner}/stack.yaml stack exec -- realigner -i {params.it} {wildcards.ref} < {params.dir}/{wildcards.ref}/hoge | samtools view -Sb | samtools sort >  {output}
#         rm {params.dir}/{wildcards.ref}/hoge
#         """


rule dump_consensus:
    input:
        config['output']+"/{ref}/top_vs_reads_sorted_top"+config['top']+"_{ref}.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_realigned_count_top"+config['top']+"_it"+config['iteration']+".tsv"
    params:
        dump_dir = config['dump_consensus']
    shell:
        """
        samtools view -h {input} | STACK_YAML={params.dump_dir}/stack.yaml stack exec -- dump-consensus {wildcards.ref} -c > {output}
        """

rule make_consensus:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_realigned_count_top"+config['top']+"_it"+config['iteration']+".tsv"
    output:
        config['output']+"/{ref}/consensus_{ref}_vs_reads_realigned_top"+config['top']+"_it"+config['iteration']+".fa"
    params:
        dir = config['output'],
    shell:
        """
        python {src_dir}/count_to_consensus.py -s {wildcards.ref} -f {input} -g True > {output}
        """

rule coverage:
    input:
        config['output']+"/{ref}/top_vs_reads_sorted_top"+config['top']+"_{ref}.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_depth_top"+config['top']+"_it"+config['iteration']+".txt"
    shell:
        """
        python {src_dir}/bam-alignment_coverage.py {input} > {output}
        """

rule region:
     input:
         config['output']+"/{ref}/{ref}_vs_reads_depth_top"+config['top']+"_it"+config['iteration']+".txt"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_region_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+".txt"
     params:
         cov = config['coverage'],
     shell:
         """
         python {src_dir}/depth-high_coverage.py {input}  {params.cov} > {output}
         """

rule terminal:
    input:
        config['output']+"/{ref}/top_vs_reads_sorted_top"+config['top']+"_{ref}.bam"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_terminal_top"+config['top']+"_it"+config['iteration']+".tsv"
    shell:
        """
        python {src_dir}/bam-alignment_terminal.py {input} > {output}
        """

rule peak:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_terminal_top"+config['top']+"_it"+config['iteration']+".tsv"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_peak_top"+config['top']+"_it"+config['iteration']+"_int"+config['interval']+"_pe"+config['peak']+".txt"
    params:
        interval = config['interval'],
        peak = config['peak'],
    shell:
        """
        python {src_dir}/terminal-find_peaks.py {input} {params.interval} {params.peak} > {output}
        """

rule cut:
    input:
        region = config['output']+"/{ref}/{ref}_vs_reads_region_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+".txt",
        peak = config['output']+"/{ref}/{ref}_vs_reads_peak_top"+config['top']+"_it"+config['iteration']+"_int"+config['interval']+"_pe"+config['peak']+".txt"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_cut_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".bed"
    params:
        cut = config['cut'],
    shell:
        """
        python {src_dir}/peak-cut_region.py {input.region} {input.peak} {params.cut} > {output}
        """

rule gapped_result:
     input:
         bed = config['output']+"/{ref}/{ref}_vs_reads_cut_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".bed",
         consensus = config['output']+"/{ref}/consensus_{ref}_vs_reads_realigned_top"+config['top']+"_it"+config['iteration']+".fa"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_result_gapped_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
     shell:
         """
         if [ -s {input.bed} ] ;then
            seqkit subseq --bed {input.bed} {input.consensus} > {output}
         else
            touch {output}
         fi
         """

rule each_result:
    input:
        config['output']+"/{ref}/{ref}_vs_reads_result_gapped_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
    output:
        config['output']+"/{ref}/{ref}_vs_reads_result_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
    shell:
        """
        python {src_dir}/remove_gap_from_fasta.py {input} > {output}
        """

rule repeat:
     input:
         expand(config['output']+"/{ref}/{ref}_vs_reads_result_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa",ref=top_id)
     output:
         config['output']+"/repeat_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
     shell:
         """
         touch {output}
         cat {input} >> {output}
         """

# rule rep_vs_rep:
#     input:
#         config['output']+"/repeat_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
#     output:
#         config['output']+"/rep_vs_rep_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".paf"
#     threads: 16
#     shell:
#         """
#         minialign -O paf {input} {input}  -t{threads} > {output}
#         """

# rule rep_vs_read:
#     input:
#         config['output']+"/repeat_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
#     output:
#         strong = config['output']+"/strong_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".paf",
#         weak = config['output']+"/weak_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".paf"
#     params:
#         data = config['data']
#     threads: 16
#     shell:
#         """
#         minialign -O paf {params.data} {input}  -t{threads} > {output.strong}
#         minialign -O paf {params.data} {input}  -t{threads} > {output.weak}
#         """
# rule read_mask:
#     input:
#         strong = config['output']+"/strong_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".paf",
#         repeat = config['output']+"/repeat_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
#     output:
#         config['output']+"/read_masked_top"+config['top']+"_it"+config['iteration']+"_cov"+config['coverage']+"_int"+config['interval']+"_pe"+config['peak']+"_cut"+config['cut']+".fa"
#     params:
#         data = config['data']
#     shell:
#         """
#         python {src_dir}/read_mask.py {params.data} {input.strong} {output}
#         """
