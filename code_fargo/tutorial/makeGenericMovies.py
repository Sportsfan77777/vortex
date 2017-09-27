"""
usage: makeGenericMovies.py [-h] [--dir DIRECTORY] [--type MOVIE_CHOICE]
                            [--save_name MOVIE_NAME]
                            [--range FILE_RANGE FILE_RANGE FILE_RANGE]
                            [--fps FPS]

optional arguments:
  -h, --help            show this help message and exit
  --dir DIRECTORY       location of files and output (default: current)
  --type MOVIE_CHOICE   select file prefix (choices: density (default),
                        polarDensity)
  --save_name MOVIE_NAME
                        select movie suffix (default: movie)
  --range START_FRAME END_FRAME RATE
                        range(start, end, rate) (default: [0, 100, 1])
  --fps FPS             movie frames per second (default: 5)
"""

import subprocess
import time

import argparse

### Movie Options ###
movie_dictionary = {}
movie_dictionary["density"] = "densityMap_"
movie_dictionary["polarDensity"] = "polarDensityMap_"

### Input Parameters ###

def new_argument_parser(description = "Manage movie input parameters."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('--dir', dest = "directory", default = ".",
                         help = 'location of files and output (default: current)')
    parser.add_argument('--type', dest = "movie_choice", default = "density",
                         help = 'select file prefix (choices: density (default), polarDensity)')
    parser.add_argument('--save_name', dest = "movie_name", default = "movie",
                         help = 'select movie suffix (default: movie)')

    # File Range
    parser.add_argument('--range', dest = "file_range", type = int, nargs = 3, default = [0, 100, 1],
                         help = 'range(start, end, rate) (default: [0, 100, 1])')

    # Movie Parameters
    parser.add_argument('--fps', dest = "fps", type = int, default = 5,
                         help = 'movie frames per second (default: 5)')

    return parser

## Parse Arguments ##
args = new_argument_parser().parse_args()

# Files
directory = args.directory
name = movie_dictionary[args.movie_choice]
movie_name = name + args.movie_name

# File Range
start_frame = args.file_range[0]
end_frame = args.file_range[1]
rate = args.file_range[2]

# Movie Parameters
fps = args.fps

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
    path = "%s/%s%s.png" % (directory, tmp_name, "%04d")
    output = "%s/%s.mp4" % (directory, movie_name)

    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    subprocess.Popen(command.split())

### Make Movies ###

# Parse Arguments
o, arguments = new_option_parser().parse_args()

if o.sample:
    sample = range(o.s_start, o.s_end, o.s_rate)
    for frame_i in sample:
        make_plot(frame_i)
elif o.frame is not None:
    make_plot(o.frame, show = True)
else:
    print "Must specify frame(s)!"

# Re-number to consecutive frames
base_path = directory
base_name = name

old_range = range(start_frame, end_frame + 1, rate)
new_range = range(len(old_range))

renumber(old_range, new_range, base_path, base_name)

# Make Movies!
make_movies()

# Delete temporary files
time.sleep(3) # 3 seconds --- (don't delete files before movie is created)
delete_tmp_files(new_range, base_path, base_name)

