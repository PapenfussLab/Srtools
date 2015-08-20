#!/usr/bin/env python

"""
fastq_split.py [-n|--num_files N_FILES] <input filename> <output directory>
"""

import os
import sys
import math
from srt.fastq import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--num_files", dest="num_files",
    help="Number of output files", type="int", default=5)
(options, args) = parser.parse_args()

input_filename = args[0]
output_directory = args[1]

if options.num_files<=0:
    print "Number of files must be > 0"
    sys.exit(-1)
num_places = 1+int(math.log10(options.num_files))

# Define output filename format
_ = os.path.split(input_filename)[-1]
base = os.path.splitext(_)[0]
ext = os.path.splitext(_)[1]
format = os.path.join(output_directory, "%s_%%0.%ii%s") % (base, num_places, ext)

# Open files
output_file = []
for i in xrange(options.num_files):
    output_filename = format % (i+1)
    output_file.append(FastqFile(output_filename, "w"))

# Split reads
i = 0
for h,s,q in FastqFile(input_filename):
    output_file[i].write(h,s,q)
    i = (i+1) % options.num_files

# Close files
for _ in output_file:
    _.close()
