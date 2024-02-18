"""
usage: makeGenericMovies.py [-h] [--dir DIRECTORY]
                            [--type {polarDensity,azimuthalDensity,density}]
                            [--id NAME_ID] [-v VERSION]
                            [--save_dir SAVE_DIRECTORY]
                            [--save_name MOVIE_NAME] [--fps FPS]
                            frames frames frames

positional arguments:
  frames                select frame range(start, end, rate)

optional arguments:
  -h, --help            show this help message and exit
  --dir DIRECTORY       location of files and output (default: current)
  --type {polarDensity,azimuthalDensity,density}
                        select file prefix (default: density)
  --id NAME_ID          id number for files -- useful for varying plot
                        parameters (default: None)
  -v VERSION            version number for files -- useful for varying plot
                        parameters (default: None)
  --save_dir SAVE_DIRECTORY
                        location of files and output (default: current)
  --movie MOVIE_NAME
                        select movie suffix (default: movie)
  --fps FPS             movie frames per second (default: 5)
"""

import os, shutil, subprocess
import argparse

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Manage movie input parameters."):
    parser = argparse.ArgumentParser()

    # Frames
    parser.add_argument('frames', type = int, nargs = 3,
                         help = 'select frame range(start, end, rate)')

    # Files
    parser.add_argument('--dir', dest = "directory", default = ".",
                         help = 'location of files and output (default: current)')
    parser.add_argument('--stem', dest = "movie_choice", default = "densityMap_",
                         help = 'select file prefix (default: density)')
    parser.add_argument('--id', dest = "name_id", type = int, default = None,
                         help = 'id number for files -- useful for varying plot parameters (default: None)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number for files -- useful for varying plot parameters (default: None)')

    parser.add_argument('--save_dir', dest = "save_directory", default = None,
                         help = 'location of files and output (default: file directory)')
    parser.add_argument('--movie', dest = "movie_name", default = None,
                         help = 'select movie suffix (default: cwd-movie)')

    # Movie Parameters
    parser.add_argument('--fps', dest = "fps", type = int, default = 5,
                         help = 'movie frames per second (default: 5)')
    parser.add_argument('--bit', dest = "bit_rate", type = int, default = 1000,
                         help = 'bit rate in kbits/s (default: 1000)')

    # Speed up on El Gato
    parser.add_argument('--fast_off', dest = "speed_up", action = 'store_false', default = True,
                         help = 'speed up on elgato w/ x264 architecture (default: do it!)')

    # Command for ffmpeg
    parser.add_argument('-t', dest = "tw2", action = 'store_true', default = False,
                         help = 'call correct ffmpeg command on tw2 (default: do not do it!)')

    return parser

## Parse Arguments ##
args = new_argument_parser().parse_args()

# Frames
start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
frame_range = range(start, end + 1, rate)

# Files
directory = args.directory
if args.save_directory is not None:
   save_directory = args.save_directory
else:
   save_directory = directory

name = args.movie_choice
if args.name_id is not None:
   name = "id%04d_%s" % (args.name_id, name)

if args.version is not None:
   name = "v%04d_%s" % (args.version, name)
if args.movie_name == None:
   cwd = os.getcwd().split("/")[-2]
   args.movie_name = "%s-movie" % cwd
movie_name = name + args.movie_name

# Movie Parameters
fps = args.fps
bit_rate = args.bit_rate

# Speed Parameter
speed_up = args.speed_up

# Command
ffmpeg = "ffmpeg"
if args.tw2:
   ffmpeg = "/home/mhammer44444/ffmpeg/bin/ffmpeg"

###############################################################################

### Helper Functions ###

def renumber(old_range, new_range, base_path, base_name):
    """ re-number files to be consecutive with 'tmp_' prefix """
    for (old_number, new_number) in zip(old_range, new_range):
        old_file = "%s/%s%04d.png" % (base_path, base_name, old_number)
        tmp_file = "%s/tmp_%s%04d.png" % (base_path, base_name, new_number)
        shutil.copyfile(old_file, tmp_file)

def delete_tmp_files(tmp_range, base_path, base_name):
    """ delete all 'tmp' files created with renumber """
    for number in tmp_range:
        tmp_file = "%s/tmp_%s%04d.png" % (base_path, base_name, number)
        os.remove(tmp_file)

def make_movies():
    """ terminal ffmpeg command """
    ### Make Movie Command (Input: frames per second; path; output_name) ###
    path = "%s/tmp_%s%s.png" % (directory, base_name, "%04d")
    output = "%s/%s.mp4" % (save_directory, movie_name)

    command = "%s -f image2 -r %d -i %s -b:v %dk -vcodec mpeg4 -y %s" % (ffmpeg, fps, path, bit_rate, output)
    if speed_up:
       command += ' -c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"'
    process = subprocess.Popen(command.split())
    process.wait()

###############################################################################

### Make Movies ###

# Re-number to consecutive frames
base_path = directory; base_name = name
new_range = range(len(frame_range))

renumber(frame_range, new_range, base_path, base_name)

# Make Movies!
make_movies()

# Delete temporary files
delete_tmp_files(new_range, base_path, base_name)

