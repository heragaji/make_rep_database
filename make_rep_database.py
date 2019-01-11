import sys
import click
import subprocess
import json
import collections as cl
import os

@click.command()
@click.option('-d','--dir', help='outout directory')
@click.option('-f','--file', help='input file')
@click.option('-r','--realigner', help='realigner directory')
@click.option('-s','--dump_consensus', help='dump-consensus directory')
@click.option('-t','--top', help='the number of top reads')
@click.option('-i','--iteration', help='the maximum number of iteration of realigner')
@click.option('-c','--coverage', help='threshold of depth')
@click.option('-k','--interval', help='threshold of peak interval')
@click.option('-p','--peak', help='threshold of peak coverage')
@click.option('-l','--cut', help='threshold of units length')

def command(dir,file,realigner,dump_consensus,top,iteration,coverage,interval,peak,cut):
    os.chdir(dir)
    subprocess.call("seqkit sort -lr {reads_data} | seqkit head -n {top} > top.fa".format(reads_data=file,top=top),shell = True)
    top_id = subprocess.check_output('seqkit seq -ni top.fa',shell = True).decode(encoding='utf-8').split()
    for id in top_id:
        os.makedirs(dir+"/"+id, exist_ok=True)
    config = cl.OrderedDict()
    config["output"] = dir
    config["data"] = file
    config["realigner"] = realigner
    config["dump_consensus"] = dump_consensus
    config["top"] = top
    config["iteration"] = iteration
    config["coverage"] = coverage
    config["interval"] = interval
    config["peak"] = peak
    config["cut"] = cut
    config["ref"] = top_id
    config["src_dir"] = os.path.dirname(__file__)+"/src"
    f = open('config.json','w')
    json.dump(config,f,indent=4)
    f.close()
    os.chdir(os.path.dirname(__file__))
    sn_log = subprocess.check_output('snakemake --directory '+ dir,shell = True).decode(encoding='utf-8')
    log = open(dir+'/make_rep_database.log','w')
    log.write(sn_log)
    log.close()


if __name__ == "__main__":
    command()

