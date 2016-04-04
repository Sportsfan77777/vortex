"""
make movies
must supply path, everything, etc.

makes use of renumber_to_tmp.py
"""

import numpy as np
import subprocess
import time

from renumber_to_tmp import tmp_prefix
from renumber_to_tmp import renumber
from renumber_to_tmp import delete_tmp_files

### Movie Commands ###
def make_movies():
    # Movie Parameter
    fps = 1

    path = "%s/%s%s.png" % (save_directory, tmp_name, "%03d")
    output = "%s/%s.mov" % (save_directory, movie_name)

    ### Make Movie Command ###
    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)

#### Paths + Names ###
save_directory = "."
name = "zoom_densityMap_"
tmp_name = "%s%s" % (tmp_prefix, name)
#movie_name = name
movie_name = name + "upTo1400"

# Re-number if necessary
base_path = save_directory
base_name = name

old_range = np.linspace(0, 1400, 15)
new_range = range(15)
renumber(old_range, new_range, base_path = base_path, base_name = base_name)

# Make Movies
make_movies()

# Delete files if necessary
time.sleep(3) # 2 seconds (don't delete files before movie is created)
delete_tmp_files(new_range, base_path = base_path, base_name = base_name)