import os

configfile:
    "config.json"


rule all:
    input:
        config['output']+"/repeat.fa"

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
         dir = config['output']
     shell:
         """
         python {src_dir}/count_to_consensus.py -s {wildcards.ref} -f {input} -g True > {output}
         """

rule coverage:
     input:
         config['output']+"/{ref}/{ref}_vs_reads_realigned.bam"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_depth.txt"
     shell:
         """
         python {src_dir}/bam-alignment_coverage.py {input} > {output}
         """

rule region:
     input:
         config['output']+"/{ref}/{ref}_vs_reads_depth.txt"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_region.txt"
     params:
         cov = config['coverage']
     shell:
         """
         python {src_dir}/depth-high_coverage.py {input}  {params.cov} > {output}
         """

rule terminal:
     input:
         config['output']+"/{ref}/{ref}_vs_reads_realigned.bam"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_terminal.tsv"
     shell:
         """
         python {src_dir}/bam-alignment_terminal.py {input} > {output}
         """

rule peak:
     input:
         config['output']+"/{ref}/{ref}_vs_reads_terminal.tsv"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_peak.txt"
     params:
         interval = config['interval'],
         peak = config['peak']
     shell:
         """
         python {src_dir}/terminal-find_peaks.py {input} {params.interval} {params.peak} > {output}
         """

rule cut:
     input:
         region = config['output']+"/{ref}/{ref}_vs_reads_region.txt",
         peak = config['output']+"/{ref}/{ref}_vs_reads_peak.txt"
     output:
         config['output']+"/{ref}/{ref}_vs_reads_cut.bed"
     params:
         cut = config['cut']
     shell:
         """
         python {src_dir}/peak-cut_region.py {input.region} {input.peak} {params.cut} > {output}
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
     shell:
         """
         python {src_dir}/remove_gap_from_fasta.py {input} > {output}
         """

rule repeat:
     input:
         expand(config['output']+"/{ref}/{ref}_vs_reads_result_gapped.fa",ref=config['ref'])
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
     shell:
         """
         python {src_dir}/cal_cluster_size.py {input} > {output}
         """
