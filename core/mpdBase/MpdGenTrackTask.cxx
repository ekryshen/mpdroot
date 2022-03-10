#include "MpdGenTrackTask.h"

MpdGenTrackTask::MpdGenTrackTask() :
fTracksInfo(nullptr) {
    fTracksInfo = new TClonesArray("MpdGenTrack");
}

MpdGenTrackTask::~MpdGenTrackTask() {
    delete fTracksInfo;
}

InitStatus MpdGenTrackTask::Init() {
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->Register("GenTracks", "GenTracks_", fTracksInfo, kTRUE);

    return kSUCCESS;
}




