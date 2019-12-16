def findElectrons(opts):
    
    ### Load Python modules

    import os
    import math
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
    h_energyCut_noTrack = TH1F("h_energyCut_noTrack","all particle energy - 20 GeV cut (NO TRACK)",nBins,eBinning)
    h_energyCut_Track = TH1F("h_energyCut_Track","all particle energy - 20 GeV cut (TRACK)",nBins,eBinning)
    
    ##BGO
    h_energyBGOl=[]  #energy of BGO vertical layer (single vertical plane)
    for BGO_idxl in range(14):
        histoName = "h_energyBGOl_" + str(BGO_idxl)
        histoTitle = "BGO energy deposit layer " + str(BGO_idxl)
        tmpHisto = TH1F(histoName,histoTitle,1000,0,1e+6)
        h_energyBGOl.append(tmpHisto)

    h_energyBGOb = [] #energy of BGO lateral layer (single bars of a plane)
    h_BGOb_maxEnergyFraction = [] #fraction of the maximum released energy for each bar on each layer of the BGO calorimeter

    for BGO_idxl in range(14):
        tmp_eLayer = []
        for BGO_idxb in range(23):
            histoName = "h_energyBGOl_" + str(BGO_idxl) + "_BGOb_" + str(BGO_idxb)
            histoTitle = "BGO energy deposit layer " + str(BGO_idxl) + " bar " + str(BGO_idxb)
            tmpHisto = TH1F(histoName,histoTitle,1000,0,1e+6)
            tmp_eLayer.append(tmpHisto)
            
        maxhistoName = "h_BGO_maxEnergyFraction_l_" + str(BGO_idxl)
        maxhistoTitle = "fraction of the maximum released energy layer " + str(BGO_idxl)
        tmpMaxHisto = TH1F(maxhistoName,maxhistoTitle,100,0,1)
        h_BGOb_maxEnergyFraction.append(tmpMaxHisto)
        h_energyBGOb.append(tmp_eLayer)

    h_BGOl_maxEnergyFraction = TH1F("h_BGOl_maxEnergyFraction","Fraction of the maximum released energy",100,0,1)

    h_thetaBGO = TH1F("h_thetaBGO","theta BGO",100,0,90)

    ##STK

    h_STK_nTracks = TH1F("h_STK_nTracks","number of tracks",1000,0,1000)
    h_STK_trackChi2norm = TH1F("h_STK_trackChi2norm","\chi^2/n track",100,0,200)
    h_STK_nTracksChi2Cut = TH1F("h_STK_nTracksChi2Cut","number of tracks (\chi^2 cut)",1000,0,1000)
        
    h_stk_cluster_XvsY = []
    for iLayer in range(6):
        hName = 'h_stkCluster_XvsY_l_'+str(iLayer)
        hTitle = 'cluster X vs Y - plane '+str(iLayer)
        tmpHisto = TH2F(hName,hTitle,1000,-500,500,1000,-500,500)
        h_stk_cluster_XvsY.append(tmpHisto)

    h_ThetaSTK = TH1F("h_ThetaSTK","theta STK",100,0,90)
    h_deltaTheta = TH1F("h_deltaTheta","\Delta theta",500,-100,100)
    
    h_resX_STK_BGO = TH1F("h_resX_STK_BGO","BGO/STK residue layer X",200,-1000,1000)
    h_resY_STK_BGO = TH1F("h_resY_STK_BGO","BGO/STK residue layer Y",200,-1000,1000)

    h_imapctPointSTK = TH2F("h_imapctPointSTK","STK impact point",1000,-500,500,1000,-500,500)

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

        #Get BGO energy deposit for each layer (vertical BGO shower profile)
        v_bgolayer  = np.array([pev.pEvtBgoRec().GetELayer(ibgo) for ibgo in range(14)])
        
        for BGO_idxl in range(14):
            h_energyBGOl[BGO_idxl].Fill(v_bgolayer[BGO_idxl])  

        #Get BGO energy deposit for each bar (lateral BGO shower profile) of each layer

        for ilay in xrange(0,14):
            v_bgolayer_bars  = np.array([pev.pEvtBgoRec().GetEdepPos(ilay,ibar) for ibar in xrange(0,23)])
            #Fraction of the maximum energy deposit of the particle crossing the BGO on a certain layer (single bars)
            h_BGOb_maxEnergyFraction[ilay].Fill(np.max(v_bgolayer_bars)/1000./etot)
            for idx_BGOb in range (23):
                h_energyBGOb[ilay][idx_BGOb].Fill(v_bgolayer_bars[idx_BGOb])
            

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

        tgZ = math.atan(np.sqrt( (pev.pEvtBgoRec().GetSlopeXZ()*pev.pEvtBgoRec().GetSlopeXZ()) + (pev.pEvtBgoRec().GetSlopeYZ()*pev.pEvtBgoRec().GetSlopeYZ()) ) );
        theta_bgo = tgZ*180./math.pi

        h_thetaBGO.Fill(theta_bgo)

        #Tracks
        ntracks = pev.NStkKalmanTrack()

        if ntracks < 0:
            print "\nTRACK ERROR: number of tracks < 0 - ABORTING\n"
            break
        if ntracks == 0:
            h_energyCut_noTrack.Fill(etot)
        
        h_STK_nTracks.Fill(ntracks)
        h_energyCut_Track.Fill(etot)

        for iTrack in range(ntracks):
            tmpTrack = pev.pStkKalmanTrack(iTrack)
            chi2_norm = tmpTrack.getChi2()/(tmpTrack.getNhitX()+tmpTrack.getNhitY()-4)
            h_STK_trackChi2norm.Fill(chi2_norm)

            if chi2_norm > 25: 
                continue
        
            h_STK_nTracksChi2Cut.Fill(ntracks)

            l0ClusterX = l0ClusterY = False

            for iCluster in range(tmpTrack.GetNPoints()):
                clux = tmpTrack.pClusterX(iCluster)
                cluy = tmpTrack.pClusterY(iCluster)
                if clux and clux.getPlane() == 0:
                    l0ClusterX = True
                if cluy and cluy.getPlane() == 0:
                    l0ClusterY = True

                # check plot for the dead region of STK
                if(clux and cluy):
                    h_stk_cluster_XvsY[clux.getPlane()].Fill(clux.GetX(),cluy.GetY())


            if l0ClusterX == False and l0ClusterY == False:
                continue

            theta_stk =math.acos(tmpTrack.getDirection().CosTheta())*180./math.pi;

            delta_theta_STK_BGO = theta_stk - theta_bgo

            h_ThetaSTK.Fill(theta_stk)
            h_deltaTheta.Fill(delta_theta_STK_BGO)

            h_imapctPointSTK.Fill(tmpTrack.getImpactPoint().x(),tmpTrack.getImpactPoint().y())
                
            h_resX_STK_BGO.Fill(tmpTrack.getImpactPoint().x() - (pev.pEvtBgoRec().GetInterceptXZ() + tmpTrack.getImpactPoint().z() * pev.pEvtBgoRec().GetSlopeXZ()))
            h_resY_STK_BGO.Fill(tmpTrack.getImpactPoint().y() - (pev.pEvtBgoRec().GetInterceptYZ() + tmpTrack.getImpactPoint().z() * pev.pEvtBgoRec().GetSlopeYZ()))
    
            if abs(theta_stk - theta_bgo) > 25:
                continue
                    
            #Now I need to select the track to extract the charge from the PSD and STK

            #Up to now I have the theta extracted from the BGO and STK

                        





    ### Writing output files to file

    if opts.data:

        tf_skim = TFile(opts.outputFile,"RECREATE")

        h_energy_all.Write()
        h_energyCut.Write()
        h_energyCut_SAAcut.Write()
        h_energyCut_noTrack.Write()
        h_energyCut_Track.Write()

        for BGO_idxl in range(14):
            h_energyBGOl[BGO_idxl].Write()
            h_BGOb_maxEnergyFraction[BGO_idxl].Write()
            for BGO_idxb in range(23):
                h_energyBGOb[BGO_idxl][BGO_idxb].Write()
        
        h_thetaBGO.Write()
        h_BGOl_maxEnergyFraction.Write()
        h_terrestrial_lat_vs_long.Write()

        h_STK_nTracks.Write()
        h_STK_trackChi2norm.Write()
        h_STK_nTracksChi2Cut.Write()

        for iLayer in range(6):
            h_stk_cluster_XvsY[iLayer].Write()

        h_ThetaSTK.Write()
        h_deltaTheta.Write()

        h_imapctPointSTK.Write()
        h_resX_STK_BGO.Write()
        h_resY_STK_BGO.Write()

        tf_skim.Close()

    
        