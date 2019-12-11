#!/usr/bin/env python

import os
import sys
from datetime import datetime
from argparse import ArgumentParser

start=datetime.now()

def main(args=None):

    ### Parsing options
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="electron skimmer")

    parser.add_argument("-i","--input",dest='input', help='single input file to process ')
    parser.add_argument("-l","--list",dest='list', help='list of files to process ')
    parser.add_argument("-v","--verbose", action='store_true', default=False, dest='verbose', help='run in high verbosity mode')
    parser.add_argument("-mc",default=False, action="store_true", dest='mc', help='use this flag for MC data')
    #parser.add_argument("-q","--quiet", action='store_true', default=False, dest='quiet', help='suppress a lot of output, quiet mode')
    parser.add_argument("-o","--outputDir",default="skimmed_data/", type=str, dest='outputDir', help='name of output directory')
    parser.add_argument("-of","--outputFile",default="skimmed_data/skimmed.root", type=str, dest='outputFile', help='name of output root file')
    parser.add_argument("-d","--debug", action='store_true', default=False, dest='debug', help='activate debig mode to select a small subsample of events for the analysis (1000)')
    parser.add_argument("-data","--data",default=True, action="store_true", dest='data', help='save analysis results in a root file')

    opts = parser.parse_args(args)
    
    #Setting output diretcory
    if not os.path.isdir(opts.outputDir):
        try:
            os.mkdir(opts.outputDir)
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