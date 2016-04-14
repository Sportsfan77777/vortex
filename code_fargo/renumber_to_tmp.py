"""
re-number to tmp files
(useful for movies)

Usage:
import renumber, delete_tmp_files

"""

import os
import shutil

tmp_prefix = "tmp_"

def renumber(old_range, new_range, base_path = None, base_name = None):
    """ cp old files to new_files with 'tmp' prefix """
    for (old_number, new_number) in zip(old_range, new_range):
        old_file = "%s/%s%04d.png" % (base_path, base_name, old_number)
        tmp_file = "%s/%s%s%04d.png" % (base_path, tmp_prefix, base_name, new_number)
        shutil.copyfile(old_file, tmp_file)

def delete_tmp_files(tmp_range, base_path = None, base_name = None):
    """ delete all 'tmp' files in a particular range """
    for number in tmp_range:
        tmp_file = "%s/%s%s%04d.png" % (base_path, tmp_prefix, base_name, number)
        os.remove(tmp_file)
