def findElectrons(opts):
    
    ### Load Python modules

    import os
    import numpy as np
    from array import array
    from os.path import isdir, abspath

    ### Load ROOT modules
    from ROOT import TClonesArray, TFile, TTree, gSystem, gROOT, AddressOf
    from ROOT import TH2F, TH1F, TMath, TGraphAsymmErrors

    ###Load DAMPE libs

    gSystem.Load("libDmpEvent.so")
    gSystem.Load("libDmpEventFilter.so")
    
    gSystem.Load("libDmpKernel.so")
    gSystem.Load("libDmpService.so")

    ###Load DAMPE modules

    from ROOT import DmpChain, DmpEvent, DmpFilterOrbit, DmpPsdBase, DmpCore
    from ROOT import DmpSvcPsdEposCor, DmpVSvc   #DmpRecPsdManager
    import DMPSW

    gROOT.SetBatch(True)

    ############################# Searching for electrons

    ####### Reading input files

    #Creating DAMPE chain for input files
    dmpch = DmpChain("CollectionTree")
    
    #Reading input files
    if not opts.input:
        files = [f.replace("\n","") for f in open(opts.list,'r').readlines()]
        for ifile, f in enumerate(files):
            DMPSW.IOSvc.Set("InData/Read" if ifile == 0 else "InData/ReadMore",f)
            if os.path.isfile(f):
                dmpch.Add(f)
                if opts.verbose:
                    print('\nInput file read: {} -> {}'.format(ifile,f))
    else:
        DMPSW.IOSvc.Set("InData/Read",opts.input)
        if os.path.isfile(opts.input):
            dmpch.Add(opts.input)
            if opts.verbose:
                print('\nInput file read: {}'.format(opts.input))
    
    #Defining the total number of events
    nevents = dmpch.GetEntries()

    if opts.verbose:
        print('\nTotal number of events: {}'.format(nevents))
        print("\nPrinting the chain...\n")
        dmpch.Print()
    
    ####### Setting the output directory to the chain
    dmpch.SetOutputDir(abspath(opts.outputDir),"electrons")

    ####### Processing input files

    ###Histos

    #Defining log binning

    #np.logspace binning
    nBins=1000
    eMax=6
    eMin=0
    eBinning = np.logspace(eMin, eMax, num=(nBins+1))
    
    #custom binning
    ''' 
    nBins = 1000
    eMin=0.1
    eMax=1000000
    EDmax = []
    EDEdge = [] 
    EDstepX=np.log10(eMax/eMin)/nBins
    for iedge in range(0, nBins):
        EDEdge.append(eMin*pow(10,iedge*EDstepX))
        EDmax.append(eMin*pow(10,(iedge+1)*EDstepX))
    EDEdge.append(EDmax[-1])
    Edges= array('d',EDEdge) # this makes a bound array for TH1F
    '''

    h_energy = TH1F("h_energy","h_energy",nBins,eBinning)

    ###

    #Filtering for SAA
    if not opts.mc:
        DMPSW.IOSvc.Set("OutData/NoOutput", "True")
        DMPSW.IOSvc.Initialize()
        pFilter = DmpFilterOrbit("EventHeader")
        pFilter.ActiveMe()
    
    #Starting loop on files

    if opts.debug:
        if opts.verbose:
            print('\nDebug mode activated... the number of chain events is limited to 1000')
        nevents = 1000

    for iev in xrange(0,nevents):
        if opts.mc:
            DmpVSvc.gPsdECor.SetMCflag(1)
        pev=dmpch.GetDmpEvent(iev)
        etot=pev.pEvtBgoRec().GetTotalEnergy()/1000.
        h_energy.Fill(etot)
    
    if opts.data:
        tf_skim = TFile(opts.outputFile,"RECREATE")
        h_energy.Write()
        tf_skim.Close()

    
        