"""
reads a title to a directory

Usage:
(1) from readTitle import readTitle
(2) python readTitle.py <=== prints title

title is stored in 'title.p'
"""

import pickle

def readTitle():
    fn = "title.p"
    try:
        title = pickle.load(open(fn, "rb"))
        return "[%s]" % title
    except:
        return None


if __name__ == "__main__":
    title = readTitle()
    print "Title: [%s]" % title