def analysisParamters(opts):
    print("\n\n ********* Analysis parameters *********\n")
    
    if not opts.input:
        print('Input file list: {}'.format(opts.list))
    else:
        print('Input single file: {}'.format(opts.input))

    print('Output directory: {}'.format(opts.output))
    print('Verbosity: {}'.format(opts.verbose))

    print("\n ********* Analysis parameters *********\n\n")