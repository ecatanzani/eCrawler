#!/usr/bin/env python

import os
import sys
from datetime import datetime
from argparse import ArgumentParser

start=datetime.now()

def main(args=None):

    ### Parsing options
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="electron skimmer")

    parser.add_argument("-i","--input",dest='input', help='the input files to be processed ')
    parser.add_argument("-v","--verbose", action='store_true', default=False, dest='verbose', help='run in high verbosity mode')
    #parser.add_argument("-q","--quiet", action='store_true', default=False, dest='quiet', help='suppress a lot of output, quiet mode')
    parser.add_argument("-o","--output",default="skimmed_data/", type=str, dest='output', help='name of output directory')

    opts = parser.parse_args(args)
    
    #Setting output diretcory
    if not os.path.isdir(opts.output):
        try:
            os.mkdir(opts.output)
        except OSError:
            print ("Creation of the output directory %s failed" % opts.output)

    #Load analysis functions
    sys.path.append("Stuff")
    from Inspector import findElectrons
    from runParameters import analysisParamters

    if(opts.verbose):
        analysisParamters(opts)

    findElectrons(opts)

if __name__ == "__main__":
    main()
    print datetime.now()-start