#ifndef MPDEMCCLUSTER_H
#define MPDEMCCLUSTER_H 1

#include "FairTask.h"
#include "TVector3.h"

using namespace std;

namespace MpdEmcMath {

  template <typename T>
  class make_vector {
  public:
   typedef make_vector<T> my_type;
   my_type& operator<< (const T& val) {
   data_.push_back(val);
   return *this;
  }
  operator std::vector<T>() const {
    return data_;
  }
  private:
   std::vector<T> data_;
  };

  inline Float_t GetPhiDiff(Float_t phi1, Float_t phi2) {
   Float_t a = phi1 - phi2;
   a += (a>180) ? -360 : (a<-180) ? 360 : 0;
   return a;
  }

}

class MpdEmcCluster : public TObject{
public:

/** Default constructor **/
    MpdEmcCluster();

/** Constructor with hit parameters **/
    MpdEmcCluster(Float_t energy, Float_t time, TVector3 pos);

/** Destructor **/
    virtual ~MpdEmcCluster();

    void Print(const Option_t* opt = 0) const;

// Calculate claster radius

    Float_t ComputeClusterRadius(TVector3 fPos, vector<Float_t> xC, vector<Float_t> yC,
                                 vector<Float_t> zC, vector<Float_t> eH);

// Set and get cluster parameters      

    Int_t GetFlag() const {
        return fFlag;
    };

    Float_t GetE() const {
        return fE;
    };

    Float_t GetTime() const {
        return fTime;
    };

    Float_t GetPhi() const {
        return fPhi;
    };

    Float_t GetRho() const {
        return fRho;
    };

    Float_t GetZ() const {
        return fZ;
    };

    Float_t GetRad() const{
	return fRad;
    };

    Float_t GetNHits() const{
	return fNHits;
    };

    void SetFlag(Int_t flag)  {
        fFlag = flag;
    };

    void SetEnergy(Float_t energy) {
	fE = energy; 
    };

    void SetTime(Float_t time) {
	fTime = time;
    };

    void SetPhi(Float_t phi) {
        fPhi = phi;
    }; 

    void SetRho(Float_t rho) {
        fRho = rho;
    };

    void SetZ(Float_t z) {
        fZ = z;
    };

    void SetRad(Float_t rad) {
	fRad = rad;
    };

    void SetNHits(UInt_t n) {
	fNHits = n;
    };

protected:

    // base digi-part
    Int_t fFlag;

    Float_t fE; // cluster energy
    Float_t fTime; // cluster time
    Float_t fPhi;// phi - angle of cluster
    Float_t fRho;// radius of custer
    Float_t fZ;// z-coordinate of cluster 
    Float_t fRad;// cluster radius
    UInt_t fNHits;// number of hits in cluster

    ClassDef(MpdEmcCluster, 1)

};


#endif
