"""
Reads a pickle file
The pickle module can be used to store a variable of any type into a file
(This is called serialization or marshalling.)
Required Parameter: pickle filename
Example Usage:
python read_pickle.py filename.p
python read_pickle.py -i list.p
Options:
  -h, --help  show this help message and exit
  -i          read files one entry at a time (only applicable for iterables)
Polish: 6
"""

import sys
import pickle as p
import collections
from optparse import OptionParser


def read_pickle(read_as_iterable):
    fn = sys.argv[-1]
    a = p.load(open(fn, 'rb'))

    # Print to stdout
    if read_as_iterable and isinstance(a, collections.Iterable):
        # requires the file is storing an iterable (e.g. a list or an array)
        for i in a:
            print i
    else:
        print a


def new_option_parser():
    """ 
    Handles input
    Returns: OptionParser
    """
    parser = OptionParser()
    parser.add_option("-i", action="store_true", 
                      dest="read_as_iterable", default = False,
                      help="read files one entry at a time (only applicable for iterables)")
    return parser


######## MAIN ########
if __name__ == '__main__':
    parser = new_option_parser()
    options, args = parser.parse_args()

    # Check for arguments
    if len(sys.argv) == 1:
        print "Requires pickle filename as an argument" 
        print "Example Usage:\npython read_pickle.py file.p\npython read_pickle.py -i file.p"

    # Read Pickle File
    read_pickle(options.read_as_iterable)
