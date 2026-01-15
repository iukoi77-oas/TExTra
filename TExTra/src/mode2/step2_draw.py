import os
import subprocess as sp
from datetime import datetime

from util.logger import *
from util.SQuIRE_toolkit import *

strand_map = {
    "none": 0,
    "rf": 1,  # First-strand
    "fr": 2,  # Second-strand
    "r": 1,   # First-strand
    "f": 2    # Second-strand
}

def draw_func(args, replicate):
    verbosity=True
    pthreads = args.threads
    strandedness=strand_map.get(args.strand.lower(), 0)  # default: 0（unstranded）
    normlib = args.normlib

    ######### START TIMING SCRIPT ############
    if verbosity:
        startTime = datetime.now()
        print("start time is:" + str(startTime), file = sys.stderr)# Prints start time
        print(os.path.basename(__file__) + '\n', file = sys.stderr) #prints script name to std err
        # print("Script Arguments" + '\n' + "=================", file = sys.stderr)
        args_dict = vars(args)
        for option,arg in args_dict.items():
            print(str(option) + "=" + str(arg), file = sys.stderr) #prints all arguments to std err
        print("\n", file = sys.stderr)
    
    draw_dir = os.path.join(args.out_dir, "bigwig")
    os.makedirs(draw_dir, exist_ok=True)

    align_dir = os.path.join(args.mode1, f"alignment/{replicate}")
    infile = find_file(align_dir,".bam", replicate, 1, True)
    if not replicate:
        replicate = get_basename(infile)

    if verbosity:
        print("Making unique and total bedgraphs "+ str(datetime.now())  + "\n",file = sys.stderr)
    
    chrominfo = args.genome.replace(".fa",".chrom.sizes")
    sp.run(
            ["samtools", "faidx", "-@", str(pthreads), args.genome],
            check=True,  # 如果命令失败则抛出异常
            stderr=sp.PIPE  # 捕获错误输出
        )
    with open(chrominfo, "w") as outfile:
        sp.run(
            ["cut", "-f1,2", f"{args.genome}.fai"],  # 等价于 cut -f1,2 hg19.fa.fai
            stdout=outfile,  # 将输出重定向到文件
            check=True,
            stderr=sp.PIPE
        )

    bedgraph_list=[]
    bedgraph(infile, strandedness, draw_dir, replicate, normlib, pthreads, bedgraph_list)
    if verbosity:
        print("Making unique and total bigwigs "+ str(datetime.now())  + "\n",file = sys.stderr)    
    make_bigwig(chrominfo, bedgraph_list)
    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing outputs at "+ str(datetime.now()) + "\n",file = sys.stderr)

        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
        print('it took: ' + str(endTime-startTime) + "\n", file = sys.stderr)