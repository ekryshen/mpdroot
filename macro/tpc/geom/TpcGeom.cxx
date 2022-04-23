
#include "TpcGeom.h"

ClassImp(TpcGeom);

TpcGeom::TpcGeom() {}

TpcGeom::TpcGeom(ofstream *f)
{
   fGeoFile = f;
}

TpcGeom::~TpcGeom()
{
   delete rot;
}

//_____________________________________________________________
void TpcGeom::BuildTPC()
{

   // std::cout << points.str() << endl;

   Mpdshape *tube = new Mpdshape(fGeoFile, "tpcChamber1", "cave", "TUBE", "air");
   tube->Fill_TUBE(tpc_z, TpcInnRad, TpcOutRad);

   tube->SetPosition(0., 0., 0.);
   tube->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
   tube->SetSegment(0);
   tube->DumpToFile();

   build_wall("In");
   build_wall("Out");
   build_membrane();

   BuildSensVolume();

   BuildEC();

   delete tube;
}

//_______________________________________________________________
void TpcGeom::BuildSensVolume()
{

   //     TPC sensitive volume

   // make half of sensitive volume
   Double_t zWidth = .5 * (tpc_chamber_z - mem_thick);

   Mpdshape *tpchalf = new Mpdshape(fGeoFile, "tpc01mod", "tpcChamber1", "TUBE", "air");
   tpchalf->Fill_TUBE(zWidth, TpcInnRad + 71.4, //    check
                      TpcOutRad - 75.);
   tpchalf->SetSegment(1);
   tpchalf->SetPosition(0., 0., tpc_chamber_z - zWidth);
   tpchalf->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0.  1.");
   tpchalf->DumpWithIncrement();

   // make field cage
   FieldCage(zWidth);

   // make sectors
   BuildSensVolume_PGON();

   //  second half sens vol
   tpchalf->SetPosition(0., 0., -tpc_chamber_z + zWidth);
   tpchalf->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. -1.");
   tpchalf->DumpToFile();
}

//__________________________________________________________________
void TpcGeom::BuildSensVolume_PGON()
{
   // sensitive volume with pad plane

   Double_t zWidth = .5 * (tpc_chamber_z - mem_thick); //

   // AZ Make the sensitive volume
   Double_t shift = 0.5 * (back_plane_Al + back_plane_G10 + pad_plane);
   zWidth -= shift;
   Mpdshape *tpcsect = new Mpdshape(f, "tpc01sv", "tpc01mod#1", "PGON", "TPCmixture");
   tpcsect->Fill_PGON(-zWidth, zWidth, sv_Y_center - sv_Y, sv_Y_center + sv_Y, 15);
   tpcsect->SetSegment(1);
   tpcsect->SetPosition(0.0, 0.0, -shift); // adjust Z-position
   tpcsect->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
   tpcsect->DumpWithIncrement();
   // AZ
   BuildPadPlane();
}

//__________________________________________________________________
void TpcGeom::padpl_layers(Double_t x_small, Double_t x_large, Double_t yWidth, TList *layers)
{

   TIter     next(layers);
   Mpdshape *obj;
   Double_t  zWidth = 0;

   while (obj = (Mpdshape *)next()) {
      //    std::cout << obj->GetVolumeName();
      if (obj->GetVolumeName() == "tpc01bpG10") {
         zWidth = .5 * back_plane_G10; //  mm
      }

      if (obj->GetVolumeName() == "tpc01bpPp") {
         zWidth = .5 * pad_plane; //
      }

      if (obj->GetVolumeName() == "tpc01bpAl") {
         zWidth = .5 * back_plane_Al; //  mm
      }

      if (zWidth == 0.) {
         cout << " Something wrong";
         continue;
      }

      obj->Fill_TRAP(x_small, x_large, yWidth, zWidth);
   }
}

//__________________________________________________________________
void TpcGeom::BuildPadPlane()
{

   // Pad plane

   Double_t center_rad = sv_Y_center;
   TList   *layers     = new TList();

   Mpdshape *layer_G10 = new Mpdshape(f, "tpc01bpG10", "tpc01mod#1", "TRAP", "G10", points.str(), position.str());
   layer_G10->SetSegment(1);
   layers->Add(layer_G10);

   Mpdshape *layer_Pad = new Mpdshape(f, "tpc01bpPp", "tpc01mod#1", "TRAP", "G10", points.str(), position.str());
   layer_Pad->SetSegment(1);
   layers->Add(layer_Pad);

   Mpdshape *layer_Al = new Mpdshape(f, "tpc01bpAl", "tpc01mod#1", "TRAP", "aluminium ", points.str(), position.str());
   layer_Al->SetSegment(1);
   layers->Add(layer_Al);

   //  layers->Print();

   Double_t zSector = tpc_chamber_z - mem_thick; // z thickness of sv
   Int_t    imodule = 1;
   Double_t zCurrent, zRot;
   Double_t kx, my, phi_cent;
   Double_t kz = .0;

   for (Int_t k = 0; k < Nsect; k++) {

      phi_cent = k * phi_step;
      kx       = center_rad * TMath::Cos(phi_cent);
      my       = center_rad * TMath::Sin(phi_cent);
      zRot     = step * k - 90.;

      //    cout << " X " << kx << " Y " << my << "Z  " << kz << endl;

      if (imodule == 1) {
         //	padpl_layers(x_small, x_large, yWidth,layers);
         padpl_layers(sv_x, sv_X, sv_Y, layers);
      }

      imodule++;

      zCurrent = .5 * zSector - 0.5 * back_plane_Al;
      layer_Al->SetPosition(kx, my, zCurrent);
      layer_Al->SetRotation(zRot, 0.0, 0.0); // AZ
      layer_Al->DumpWithIncrement();

      zCurrent = zCurrent - 0.5 * (back_plane_Al + back_plane_G10);
      layer_G10->SetPosition(kx, my, zCurrent);
      layer_G10->SetRotation(zRot, 0.0, 0.0); // AZ
      layer_G10->DumpWithIncrement();

      zCurrent = zCurrent - 0.5 * (back_plane_G10 + pad_plane);
      layer_Pad->SetPosition(kx, my, zCurrent);
      layer_Pad->SetRotation(zRot, 0.0, 0.0); // AZ
      layer_Pad->DumpWithIncrement();
   }
   cout << "==== pad planes filled ====\n";
}

//_______________________________________________________________
int TpcGeom::build_membrane()
{
   //     make TPC membrane

   Mpdshape *membrane = new Mpdshape(fGeoFile, "tpc01mb", "tpcChamber1", "TUBE", "rohacellhf71");
   membrane->Fill_TUBE(mem_thick, TpcInnRad + 71., TpcOutRad - 75.);
   membrane->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
   membrane->SetPosition(0., 0., 0.);
   membrane->SetSegment(0);
   membrane->DumpToFile();
   delete membrane;
}

//_______________________________________________________________
int TpcGeom::build_wall(TString which)
{

   TString wallname = "tpc01" + which + "Wall";
   cout << wallname.Data() << endl;

   Double_t wall_thickness, radius, co2;
   wall_thickness = 2. * (2. * (thick_al + thick_tedlar) + kevlar);

   if (which == "In") {
      co2 = CO_2_gap[1];
      wall_thickness += co2;
      radius = TpcInnRad;
   } else {
      co2 = CO_2_gap[0];
      wall_thickness += co2;
      radius = TpcOutRad - wall_thickness;
   }

   cout << which << " wall thickness = " << wall_thickness << endl;

   Mpdshape *wall = new Mpdshape(fGeoFile, wallname.Data(), "tpcChamber1", "TUBE", "air");
   wall->Fill_TUBE(tpc_chamber_z, radius, radius + wall_thickness);
   wall->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
   wall->SetPosition(0., 0., 0.);
   wall->SetSegment(0);
   wall->DumpToFile();

   Double_t R_old, R_cur;

   for (Int_t section = 1; section < 3; section++) {

      if (section == 1) {
         R_old = radius;
      } else {
         R_old = R_cur;
      }

      // add  Alum
      R_cur = R_old + thick_al;

      wall->Fill_TUBE(tpc_chamber_z, R_old, R_cur);

      // std::cout << points.str() << endl;
      TString volume = "tpc01" + which + "Sct";

      volume += section;
      TString volname = volume + "Al1";

      //    std::cout << volume << endl;

      wall->SetVolumeName(volname.Data());
      wall->SetMotherVolumeName(wallname.Data());
      wall->SetMedia("aluminium");
      wall->DumpToFile();

      // tpc wall inner tube

      Double_t R_old = R_cur;
      R_cur          = R_cur + thick_tedlar;

      wall->Fill_TUBE(tpc_chamber_z, R_old, R_cur);

      // std::cout << points.str() << endl;

      volname = volume + "Tdl1";
      wall->SetVolumeName(volname.Data());
      wall->SetMedia("tedlar");
      wall->DumpToFile();

      // tpc wall inner tube kevlar

      R_old = R_cur;
      R_cur = R_cur + kevlar;

      wall->Fill_TUBE(tpc_chamber_z, R_old, R_cur);

      // std::cout << points.str() << endl;

      volname = volume + "Kvl";
      wall->SetVolumeName(volname.Data());
      wall->SetMedia("kevlar");
      wall->DumpToFile();

      // tpc  inner wall tedlar

      R_old = R_cur;
      R_cur = R_cur + thick_tedlar;

      wall->Fill_TUBE(tpc_chamber_z, R_old, R_cur);

      volname = volume + "Tdl2";
      wall->SetVolumeName(volname.Data());
      wall->SetMedia("tedlar");
      wall->DumpToFile();

      // tpc inner wall al

      R_old = R_cur;
      R_cur = R_cur + thick_al;
      wall->Fill_TUBE(tpc_chamber_z, R_old, R_cur);

      volname = volume + "Al2";
      wall->SetVolumeName(volname.Data());
      wall->SetMedia("aluminium");
      wall->SetSegment(0);
      wall->DumpToFile();

      // tpc co2 gap inner wall

      if (section == 2) continue;

      R_old = R_cur;
      R_cur = R_cur + co2;

      wall->Fill_TUBE(tpc_chamber_z, R_old, R_cur);

      volname = volume + "CO2";
      wall->SetVolumeName(volname.Data());
      wall->SetMedia("CO2");
      wall->SetSegment(0);
      wall->DumpToFile();
   }
   delete wall;

   cout << "==== " << which << " wall filled ====\n";
}

//_______________________________________________________________
void TpcGeom::FieldCage(Double_t zWidth)
{
   //================== FieldCage===================================

   Mpdshape *cage = new Mpdshape(fGeoFile, "tpc01outfc", "tpc01mod#1", "PGON", "mylar");

   // make outer fieldcage
   cage->Fill_PGON(-zWidth, zWidth, fc_out_rmin, fc_out_rmax, 7.);
   cage->SetSegment(1);
   cage->SetPosition(0., 0., 0.);
   cage->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. 1.");
   cage->DumpWithIncrement();

   // make inner fieldcage
   cage->Fill_PGON(-zWidth, zWidth, fc_in_rmin, fc_in_rmax, 7.);
   cage->SetMotherVolumeName("tpc01mod#1");
   cage->SetVolumeName("tpc01infc");
   cage->SetSegment(1);
   cage->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. 1.");
   cage->DumpWithIncrement();

   delete cage;

   std::string name[2] = {"tpc01outpin", "tpc01inpin"};

   Double_t phi_cent, kx, my, center_rad;

   Mpdshape *pin = new Mpdshape(f, "tpc01inpin", "tpc01mod#1", "TUBE", "polypropylene");
   for (Int_t j = 0; j < 2; j++) {
      center_rad = fc_r[j];

      pin->Fill_TUBE(zWidth, fc_pin_r[j] - fc_wall_thick, fc_pin_r[j]);

      pin->SetSegment(1);
      pin->SetVolumeName(name[j]);

      for (Int_t k = 0; k < Nsect; k++) {

         phi_cent = k * phi_step + fc_phi_shift; // radians
         kx       = center_rad * TMath::Cos(phi_cent);
         my       = center_rad * TMath::Sin(phi_cent);
         //      cout << " X " << kx << " Y " << my << "Z  " << kz << endl;
         pin->SetPosition(kx, my, 0.);
         pin->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. 1.");
         pin->DumpWithIncrement();
      }
   }

   delete pin;
}

//_______________________________________________________________
void TpcGeom::BuildEC()
{

   // ===========================================================

   Double_t zWidth = .5 * (tpc_z - tpc_chamber_z);

   //  right endcap
   Mpdshape *ec = new Mpdshape(fGeoFile, "tpc01ec", "tpcChamber1", "TUBE", "air");
   ec->Fill_TUBE(zWidth, TpcInnRad, TpcOutRad);
   ec->SetSegment(1);

   Double_t z_ec = zWidth + tpc_chamber_z;
   ec->SetPosition(0., 0., z_ec);
   ec->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0.  1.");
   ec->DumpWithIncrement();

   BuildFlanch(zWidth);
   BuildRibs(zWidth);

   //  left end-cup:
   ec->SetPosition(0., 0., -z_ec);
   ec->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. -1.");
   ec->DumpToFile();
   delete ec;
}

//_______________________________________________________________
void TpcGeom::BuildFlanch(Double_t zWidth)
{

   // flanch outer

   Mpdshape *flanch = new Mpdshape(fGeoFile, "tpc01of", "tpc01ec#1", "TUBE", "aluminium");
   flanch->Fill_TUBE(.5 * flanch_thickness, outer_flanch_in_rad, TpcOutRad);
   flanch->SetPosition(0., 0., .5 * flanch_thickness - zWidth);
   flanch->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. 1.");
   flanch->DumpToFile();

   //  flanch  inner

   flanch->Fill_TUBE(.5 * flanch_thickness, TpcInnRad, inner_flanch_out_rad);

   flanch->SetVolumeName("tpc01if");
   flanch->SetPosition(0., 0., .5 * flanch_thickness - zWidth);
   flanch->SetRotation(" 1. 0. 0. 0. 1. 0. 0. 0. 1.");
   flanch->DumpToFile();
}

void TpcGeom::BuildRibs(Double_t zWidth)
{
   // ==========  Ribs  ============================

   Double_t xWidth = 0.5 * rib_width_x;                                                  // x - half-dimension of block
   Double_t yWidth = (TMath::Sqrt(TpcOutRad * TpcOutRad - xWidth * xWidth) - TpcInnRad); // y - half-dimension of block
   Double_t zRibs  = .5 * rib_width_z;                                                   // z - half-dimension of block

   yWidth -= .5 * (outer_flanch_width + inner_flanch_width);

   yWidth /= 2.;

   Double_t kz = flanch_thickness + rib_pos_z - zWidth + zRibs;

   // tpc ribs

   Mpdshape *tpcrib = new Mpdshape(f, "tpc01Rib", "tpc01ec#1", "BOX", "aluminium");
   tpcrib->Fill_BOX(xWidth, yWidth, zRibs);
   tpcrib->SetSegment(1);

   Double_t kx, my, zRot, phi_cent;

   Double_t center_rad = TpcInnRad + yWidth + .5 * inner_flanch_width;
   //  cout << "center_rad_rib: " << center_rad << "\n";

   for (Int_t k = 0; k < Nsect; k++) {

      phi_cent = k * phi_step + .5 * phi_step;
      kx       = center_rad * TMath::Cos(phi_cent);
      my       = center_rad * TMath::Sin(phi_cent);
      zRot     = step * k - 90. + .5 * step;

      tpcrib->SetPosition(kx, my, kz);
      tpcrib->SetRotation(zRot, 0.0, .0);
      tpcrib->DumpWithIncrement();
   }
}

#if 0
//_________________________________________________________
void external_sector_layers(Double_t x_small, Double_t x_large,
			    Double_t yWidth, Mpdshape *layer_G10,
			    Mpdshape *layer_Cu, Mpdshape *layer_Al) {

  // G10 plane

  Double_t zWidth = .5*zG10; //  mm
  Fill_TRAP(x_small, x_large, yWidth, zWidth);

  position << 0.0 << " " << 0.0 << " " << 0.0;

  layer_G10->SetPoints(points.str());
  
  // Cu plane


  zWidth = .5*zCu; //

  Fill_TRAP(x_small, x_large, yWidth, zWidth);

  position << 0.0 << " " << 0.0 << " " << 0.0;

  layer_Cu->SetPoints(points.str());

 // Al plane

  zWidth = .5*zAl; //  mm
  Fill_TRAP(x_small, x_large, yWidth, zWidth);

  position << 0.0 << " " << 0.0 << " " << 0.0;
  layer_Al->SetPoints(points.str());

}

//__________________________________________________________________
void external_sector() {

  // Sector volume
  TVector3 CircleCenter = TVector3(0.,0.,0.);
  TVector3 LineStart = TVector3( 0.5*rib_width_x,0.,0.);
  TVector3 LineEnd = TVector3( 0.5*rib_width_x, 500.,0.);

  TVector3 crosspoint = LineCrossesCircle(CircleCenter, outer_flanch_in_rad,
					  LineStart, LineEnd);

  Double_t yWidth, zWidth, kz;
  Double_t kx,my,phi_cent;
  //  crosspoint.Print();

  Double_t x_large = (crosspoint.Y()-.5*rib_width_x/TMath::Tan(phi_step/2.))*
    TMath::Sin(phi_step/2.);
  Double_t x_small = (inner_flanch_out_rad*TMath::Sin(phi_step/2.) -
    .5*rib_width_x)/TMath::Cos(phi_step/2.) ;

  Double_t zSector = tpc_z - tpc_chamber_z; // z thickness of trap

  yWidth = .5*(x_large - x_small)/TMath::Tan(phi_step/2.);

  zWidth = .5*zSector; // 

  kz = 0.;

  //  make sector

  Fill_TRAP(x_small, x_large, yWidth, zWidth);

  position << 0.0 << " " << 0.0 << " " << 0.0;

  Double_t center_rad = inner_flanch_out_rad + yWidth;

  Mpdshape* extsect = new Mpdshape(f, "tpc01ES", "tpc01ec#1", "TRAP", "air",
				   points.str(), position.str());
  extsect->SetSegment(1);

  Int_t imodule=1;

  Double_t zCurrent,zRot;

  Mpdshape *layer_G10 = new Mpdshape(f, "tpc01ESG10", "tpc01ES#1", "TRAP", "G10",
			     points.str(), position.str());
  layer_G10->SetSegment(1);

  Mpdshape *layer_Cu = new Mpdshape(f, "tpc01ESCu", "tpc01ES#1",
			    "TRAP", "copper", points.str(), position.str());
  layer_Cu->SetSegment(1);

  Mpdshape *layer_Al = new Mpdshape(f, "tpc01ESAl", "tpc01ES#1", "TRAP",
			    "aluminium ",points.str(), position.str());
  layer_Al->SetSegment(1);
  // ========================================

    for (Int_t k = 0; k < Nsect; k++) {

      phi_cent = k*phi_step; // radians
      kx = center_rad* TMath::Cos(phi_cent);
      my = center_rad* TMath::Sin(phi_cent);
      zRot = step*k  - 90.;
      //      cout << " X " << kx << " Y " << my << "Z  " << kz << endl;
      extsect->SetPosition(kx, my, kz);
      extsect->SetRotation(zRot, 0.0, 0.0);
      extsect->DumpWithIncrement();

      //+++++++++++++++++==================
#ifdef FEESIM
      if (imodule == 1) {
	external_sector_layers(x_small, x_large, yWidth,
			       layer_G10, layer_Cu, layer_Al);
      }

      imodule++;

      zCurrent =  -.5*zSector + 0.5*zG10;

      layer_G10->SetPosition(0., 0., zCurrent);
      layer_G10->DumpWithIncrement();

      zCurrent =  zCurrent + 0.5*(zG10 + zCu);
      layer_Cu->SetPosition(0., 0., zCurrent);
      layer_Cu->DumpWithIncrement();
      
      zCurrent =  zCurrent + 0.5*(zAl + zCu);
      layer_Al->SetPosition(0., 0.,zCurrent);
      layer_Al->DumpWithIncrement();
#endif      
    }
  cout << "==== external sectors filled ====\n";
}

#endif
