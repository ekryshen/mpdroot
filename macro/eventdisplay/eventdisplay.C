R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"
void SetTasks(MpdEventManager* fMan);

// EVENT DISPLAY macro for MPD simulated data
//
// sim_file - path to the file with MC data and detector geometry
// reco_file - file with reconstructed data for simulated events
// is_online: false - use Offline Mode (manual switching of events), default; true - use Online Mode (continious viewing of events)
void eventdisplay(const char* sim_file = "$VMCWORKDIR/macro/mpd/evetest.root", const char* reco_file = "$VMCWORKDIR/macro/mpd/mpddst.root", bool is_online = false)
{
    gDebug = 0;

    // CREATE FairRunAna
    FairRunAna* fRun = new FairRunAna();

    // check file existence with MC data and detector geometry
    TString strSimFile(sim_file);
    if (!CheckFileExist(strSimFile))
    {
        cout<<endl<<"ERROR: Simulation file with detector geometry wasn't found!"<<endl;
        return;
    }

    // set source of events to display and addtiional parameters
    FairSource* fFileSource = new FairFileSource(strSimFile);

    // set parameter file with MC data and detector geometry
    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    FairParRootFileIo* parIo1 = new FairParRootFileIo();
    parIo1->open(strSimFile);
    rtdb->setFirstInput(parIo1);
    rtdb->setOutput(parIo1);
    rtdb->saveOutput();

    // add file with reconstructed data as a friend
    TString strDstFile(reco_file);
    if (CheckFileExist(strDstFile))
        ((FairFileSource*)fFileSource)->AddFriend(strDstFile);
    else
        cout<<endl<<"Warning: File with reconstructed data wasn't found!"<<endl;

    fRun->SetSource(fFileSource);

    // Create Event Manager
    MpdEventManager* fMan = new MpdEventManager();
    fMan->isOnline = is_online;

    // set output file
    fRun->SetOutputFile("ed_out.root");

    // set tasks to draw
    SetTasks(fMan);

    //if (gGeoManager)
    //    gGeoManager->SetVisLevel(3);

    // MpdEventManager::Init(Int_t visopt = 1, Int_t vislvl = 3, Int_t maxvisnds = 10000)
    fMan->Init();

    // -----   Finish   -------------------------------------------------------
    cout<<endl<<endl<<"Event Display was initialized successfully."<<endl<<endl;
    // ------------------------------------------------------------------------
}

// set FairRunAna drawing tasks depending from data source and on/offline mode
void SetTasks(MpdEventManager* fMan)
{
    Style_t pointMarker = kFullDotSmall;
    Color_t mcPointColor = kOrange, recoPointColor = kBlack;

    // draw MC points
    MpdMCPointDraw* TpcPoint = new MpdMCPointDraw("TpcPoint", mcPointColor, kDot);
    fMan->AddTask(TpcPoint);
    MpdMCPointDraw* TofPoint = new MpdMCPointDraw("TOFPoint", mcPointColor, kDot);
    fMan->AddTask(TofPoint);
    MpdMCPointDraw* ZdcPoint = new MpdMCPointDraw("ZdcPoint", mcPointColor, kFullSquare);
    fMan->AddTask(ZdcPoint);
    MpdMCPointDraw* EmcPoint = new MpdMCPointDraw("EmcPoint", mcPointColor, kDot);
    fMan->AddTask(EmcPoint);
    /* unsupported as of 12.2022
    MpdMCPointDraw* McordPoint = new MpdMCPointDraw("McordPoint", mcPointColor, kDot);
    fMan->AddTask(McordPoint); */

    // draw MC geometry tracks
    MpdMCTracks* GeoTrack = new MpdMCTracks("GeoTracks");
    fMan->AddTask(GeoTrack);

    // draw MC tracks
    //MpdMCStack* MCTrack = new MpdMCStack("MCTrack");
    //fMan->AddTask(MCTrack);

    // DST hits
    MpdHitPointSetDraw* MpdTpcHit = new MpdHitPointSetDraw("MpdTpcHit", mcPointColor, kDot);
    fMan->AddTask(MpdTpcHit);

    MpdHitPointSetDraw* MpdTofHit = new MpdHitPointSetDraw("TOFHit", mcPointColor, kDot);
    fMan->AddTask(MpdTofHit);

    // DST hits (box view)
    //MpdHitDraw* MpdTpcHit = new MpdHitDraw("MpdTpcHit", 1);
    //fMan->AddTask(MpdTpcHit);

    // DST tracks
    MpdGlobalTrackDraw* MpdGlobalTrack = new MpdGlobalTrackDraw("GlobalTracks");
    fMan->AddTask(MpdGlobalTrack);

    // draw EMC towers
    MpdEmcTowerDraw* MpdEmcTower= new MpdEmcTowerDraw("MpdEmcTower", 0.0075 /*energy threshold, GeV*/, 1);
    fMan->AddTask(MpdEmcTower);

    // draw ZDC towers
    MpdZdcTowerDraw* MpdZdcTower= new MpdZdcTowerDraw("MpdZdcTower", 0.0001 /*energy threshold, GeV*/, kFALSE, 1);
    fMan->AddTask(MpdZdcTower);

    // save EventDisplay Screenshot
    //MpdWebScreenshots* WebScreenshots = new MpdWebScreenshots("WebScreenshots", "/var/www/html/events"); // for WEB-page
    //MpdWebScreenshots* WebScreenshots = new MpdWebScreenshots("WebScreenshots","screenshots"); // folder to save the screenshots
    //WebScreenshots->SetFormatFiles(0); // 0 -.png, 1 -.jpg, 2 -.jpg and .png
    //WebScreenshots->SetMultiFiles(0); //0 - the same file (event.png), 1 - multiple files (event_nnn.png)
    //WebScreenshots->SetPort(8016); // 8016 by default
    //fMan->AddTask(WebScreenshots);
}
