def analysisParamters(opts):
    print("\n\n ********* Analysis parameters *********\n")
    
    if not opts.input:
        print('Input file list: {}'.format(opts.list))
    else:
        print('Input single file: {}'.format(opts.input))

    print('Output directory: {}'.format(opts.outputDir))
    print('Output file: {}'.format(opts.outputFile))
    print('Verbosity: {}'.format(opts.verbose))
    print('Debug: {}'.format(opts.debug))

    if opts.ms:
        print("\n MC flag activated")
    
    print("\n ***************************************\n\n")