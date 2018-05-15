"""
makes a new job script

"""


import argparse


def new_argument_parser(description = "Make a new job script."):
    parser = argparse.ArgumentParser()

    # Basic Parameters
    parser.add_argument('-c', dest = "num_cores", type = int, default = 16,
                         help = 'number of cores (default: 16)')
    parser.add_argument('-p', dest = "ptile", type = int, default = 16,
                         help = 'number of cores needed on each computer (default: 16)')

    parser.add_argument('--err', dest = "err_name", default = "err_%I",
                         help = 'job error file name (default: err_%I)')
    parser.add_argument('--out', dest = "out_name", default = "out_%I",
                         help = 'job output file name (default: out_%I)')

    parser.add_argument('-q', dest = "queue", default = "medium",
                         help = 'queue (default: medium)')
    parser.add_argument('--name', dest = "name", default = "auto",
                         help = 'queue (default: auto)')

    # Modules

    # Job
    parser.add_argument('-j', dest = "job", default = "",
                         help = 'job command (default: empty string)')
    parser.add_argument('-o', dest = "output", default = "this.out",
                         help = 'output file (default: this.out)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

if args.ptile > args.num_cores:
    args.ptile = args.num_cores

###############################################################################

### Write File ###

with open(args.fn, 'w') as f:
    f.write("#!/bin/bash\n")

    # Basic Parameters
    f.write("#BSUB -n %d\n" % args.num_cores)
    f.write("#BSUB -e %s\n" % args.err_name)
    f.write("#BSUB -o %s\n" % args.out_name)
    f.write('#BSUB -q "%s"\n' % args.queue)
    f.write("#BSUB -u mhammer\n")
    f.write("#BSUB -J %s\n" % args.name)
    if args.gpu:
        f.write("#BSUB -R gpu\n")
    f.write('#BSUB -R "span[ptile=%d]"\n' % args.ptile)

    # Line Break
    f.write("\n")

    # Modules
    if args.python:
        f.write("module load python/2.7.3\n")
    if args.fftw:
        f.write("module load fftw/2.1.5\n")
    if args.openmpi:
        f.write("module load openmpi\n")

    # Job
    if args.mpirun:
        f.write("mpirun -np %d " % args.num_cores)
    f.write("%s " % args.job)
    f.write("> %s\n" % args.output)



