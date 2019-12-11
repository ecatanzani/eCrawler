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

    #Pointing
    h_terrestrial_lat_vs_long =  TH2F("h_terrestrial_lat_vs_long","latitude vs longitude",360,0,360,180,-90,90)

    ## Energy

    h_energy_all = TH1F("h_energy_all","all particle energy",nBins,eBinning)
    h_energyCut = TH1F("h_energyCut","all particle energy - 20 GeV cut",nBins,eBinning)
    h_energyCut_SAAcut = TH1F("h_energyCut_SAAcut","all particle energy - 20 GeV cut (no SAA)",nBins,eBinning)
    ##BGO

    h_energyBGOl=[]
    for BGO_idxl in range(14):
        histoName = "h_energyBGOl_" + str(BGO_idxl)
        histoTitle = "BGO energy deposit layer " + str(BGO_idxl)
        tmpHisto = TH1F(histoName,histoTitle,1000,0,1e+6)
        h_energyBGOl.append(tmpHisto)

    h_BGOl_maxEnergyFraction = TH1F("h_BGOl_maxEnergyFraction","Fraction of the maximum released energy",100,0,1)

    ###

    ### Analysis cuts

    eCut = 50       #Energy cut in GeV

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
    
    for iev in xrange(nevents):

        if opts.mc:
            DmpVSvc.gPsdECor.SetMCflag(1)
        pev=dmpch.GetDmpEvent(iev)

        #Get latitude and longitude
        longitude = pev.pEvtAttitude().lon_geo
        latitude = pev.pEvtAttitude().lat_geo

        #Get particle total energy
        etot=pev.pEvtBgoRec().GetTotalEnergy()/1000.
        h_energy_all.Fill(etot)
        if etot < eCut:
            continue
        h_energyCut.Fill(etot)

        #Get BGO energy deposit for each layer
        v_bgolayer  = np.array([pev.pEvtBgoRec().GetELayer(ibgo) for ibgo in range(14)])
        
        for BGO_idxl in range(14):
            h_energyBGOl[BGO_idxl].Fill(v_bgolayer[BGO_idxl])  

        #Fraction of the maximum energy deposit of the particle crossing the BGO
        h_BGOl_maxEnergyFraction.Fill(np.max(v_bgolayer)/1000./etot)

        #SAA filter
        if not opts.mc:
            inSAA = pFilter.IsInSAA(pev.pEvtHeader().GetSecond())
            #inSAA = False
            if (inSAA): 
                continue
            h_energyCut_SAAcut.Fill(etot)
            h_terrestrial_lat_vs_long.Fill(longitude,latitude)

    ### Writing output files to file

    if opts.data:

        tf_skim = TFile(opts.outputFile,"RECREATE")

        h_energy_all.Write()
        h_energyCut.Write()
        h_energyCut_SAAcut.Write()

        for BGO_idxl in range(14):
            h_energyBGOl[BGO_idxl].Write()
            
        h_BGOl_maxEnergyFraction.Write()
        h_terrestrial_lat_vs_long.Write()

        tf_skim.Close()

    
        