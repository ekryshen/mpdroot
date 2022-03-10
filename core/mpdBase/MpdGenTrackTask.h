#include <FairTask.h>
#include <TClonesArray.h>
#include <FairRootManager.h>

#ifndef MPDGENTRACKTASK_H
#define MPDGENTRACKTASK_H

using namespace std;

class MpdGenTrackTask : public FairTask {
public:
    
    MpdGenTrackTask();    
    virtual InitStatus Init();      
    virtual ~MpdGenTrackTask();
    
    TClonesArray* GetTracksInfo() {
        return fTracksInfo;
    }

private:
    TClonesArray* fTracksInfo;
    
    ClassDef(MpdGenTrackTask, 1);
};
#endif 