"""
makes a new job script

"""


import argparse


def new_argument_parser(description = "Make a new job script."):
    parser = argparse.ArgumentParser()

    

    return parser

###############################################################################


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



