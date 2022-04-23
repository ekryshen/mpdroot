/******************************************************************************
 *  Version:    09-Jun-2012   (EL)
 *****************************************************************************/
class FairGeoRotation;

#include <TVector3.h>

#include <iostream>
#include <fstream>
#include <string>

class Mpdshape : public TObject {

public:
   Mpdshape();
   Mpdshape(ofstream *file, std::string volume_name, std::string mother_volume_name, std::string shape,
            std::string media, std::string points, std::string position = "0. 0. 0.",
            std::string rotation = "1. 0. 0. 0. 1. 0. 0. 0. 1.", Int_t segment = 0, Int_t mother_segment = 0);
   Mpdshape(ofstream *file, std::string volume_name, std::string mother_volume_name, std::string shape,
            std::string media);
   ~Mpdshape();

   void SetFile(ofstream *file) { fFile = file; }
   void SetVolumeName(std::string name) { fVolumeName = name; }
   void SetMotherVolumeName(std::string name) { fMotherVolumeName = name; }
   void SetShape(std::string shape) { fShape = shape; }
   void SetMedia(std::string media) { fMedia = media; }
   void SetPoints(std::string points) { fPoints = points; }

   void SetPosition(std::string position) { fPosition = position; }
   void SetPosition(Double_t x, Double_t y, Double_t z, Int_t precision = 6);
   void SetPosition(TVector3 vec);

   void SetRotation(std::string rotation) { fRotation = rotation; }
   void SetRotation(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5, Double_t p6, Double_t p7,
                    Double_t p8, Double_t p9);
   void SetRotation(Double_t z, Double_t y1, Double_t z1);

   void SetPositionRotation(std::string position, std::string rotation)
   {
      fPosition = position;
      fRotation = rotation;
   }
   void SetSegment(Int_t segment) { fSegment = segment; }
   void SetMotherSegment(Int_t segment) { fMotherSegment = segment; }

   void DumpToFile();
   void DumpWithIncrement()
   {
      DumpToFile();
      fSegment++;
   }
   static std::string itoa(Double_t x);

   std::string GetVolumeName() const { return fVolumeName; }
   Int_t       GetSegment() const { return fSegment; }
   void        Clear();
   void        Fill_TUBE(Double_t length_z, Double_t r_max, Double_t r_min);
   void        Fill_TUBS(Double_t zmin, Double_t zmax, Double_t rmin, Double_t rmax, Double_t dr_st, Double_t dr_end);
   void Fill_TRAP(Double_t dx11, Double_t dx12, Double_t dy1, Double_t dx21, Double_t dx22, Double_t dy2, Double_t dz);
   void Fill_TRAP(Double_t x, Double_t X, Double_t x_small_f, Double_t x_large_f, Double_t yW, Double_t yW2,
                  Double_t zW, Double_t zW2, Int_t ra);
   // void Fill_RECTRAP(Double_t x_small, Double_t x_large,Double_t xx_small, Double_t xx_large,
   //	 Double_t yWidth, Double_t zWidth);
   void Fill_SPHE(Double_t rzmin, Double_t rmax, Double_t thetamin, Double_t thetamax, Double_t phimin,
                  Double_t phimax);
   void Fill_PGON(Double_t zmin, Double_t zmax, Double_t rmin, Double_t rmax, Double_t phi1);
   // void Fill_BOX(Double_t x, Double_t y, Double_t z);
   void Fill_BRIK(Double_t xWidth1, Double_t xWidth2, Double_t yWidth, Double_t zWidth1, Double_t zWidth2, Int_t ra);

protected:
   std::ofstream   *fFile;
   std::string      fVolumeName;
   std::string      fMotherVolumeName;
   std::string      fShape;
   std::string      fMedia;
   std::string      fPoints;
   std::string      fPosition;
   std::string      fRotation;
   Int_t            fSegment;
   Int_t            fMotherSegment;
   FairGeoRotation *rot;

   ClassDef(Mpdshape, 0);
};
