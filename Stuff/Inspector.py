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
    files = [f.replace("\n","") for f in open(opts.input,'r').readlines()]

    for ifile, f in enumerate(files):
        DMPSW.IOSvc.Set("InData/Read" if ifile == 0 else "InData/ReadMore",f)
        print f
        if os.path.isfile(f):
            dmpch.Add(f)
            if opts.verbose:
                print ifile , f
        