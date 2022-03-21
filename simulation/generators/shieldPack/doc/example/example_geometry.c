// Example geometry for testing TShield as autonomic system.
void example_geometry() {
    double cm = 1;
    double degree = 1;
    cout << "Root, TShield example" << endl;
    gSystem->Load("libGeom.so");
    new TGeoManager("World", "Geometry for example.c root file.");

    TGeoElementTable *table = gGeoManager->GetElementTable();
    TGeoMaterial *he = new TGeoMaterial("He", table->GetElement(2), 1);
    TGeoMaterial *al = new TGeoMaterial("Al", table->GetElement(13), 1);
    TGeoMaterial *cl = new TGeoMaterial("Cl", table->GetElement(17), 1);
    TGeoMaterial *fe = new TGeoMaterial("Fe", table->GetElement(26), 1);
    TGeoMaterial *h  = new TGeoMaterial("H" , table->GetElement(1), 1);
    TGeoMaterial *pb = new TGeoMaterial("Pb", table->GetElement(82), 1);
    TGeoMaterial *vac = new TGeoMaterial("Vacuum", 0, 0, 0);

    TGeoMedium *mvac = new TGeoMedium("Vacuum", 1, vac);
    TGeoMedium *mal = new TGeoMedium("Al", 1, al);
    TGeoMedium *mhe = new TGeoMedium("He", 1, he);
    TGeoMedium *mcl = new TGeoMedium("Ar", 1, cl);
    TGeoMedium *mfe = new TGeoMedium("Fe", 1, fe);
    TGeoMedium *mh  = new TGeoMedium("H" , 1, h);
    TGeoMedium *mpb = new TGeoMedium("Pb", 1, pb);
    TGeoMixture *mhehar = new TGeoMixture("He+H+Ar", 1, 0.5);
    mhehar->AddElement(he, 0.5);
    mhehar->AddElement(h, 0.1);
    mhehar->AddElement(cl, 0.4);

    TGeoVolume *top = gGeoManager->MakeBox("Top", mfe, 100 * cm, 50 * cm, 50 * cm);

    TGeoPgon *pcon = new TGeoPcon(0, 360, 5);
    {
        pcon->DefineSection(0, 0, 0, 40 * cm);
        pcon->DefineSection(1, 20 * cm, 0, 40 * cm);
        pcon->DefineSection(2, 30 * cm, 0, 30 * cm);
        pcon->DefineSection(3, 40 * cm, 0, 40 * cm);
        pcon->DefineSection(4, 100 * cm, 0, 40 * cm);
    }
    TGeoVolume *vpcon = new TGeoVolume("PCon", pcon, mal);
    TGeoSphere *sph = new TGeoSphere(0, 30 * cm, 0, 180 * degree, 0, 360 * degree);
    TGeoVolume *vsph = new TGeoVolume("Sph", sph, mcl);
    TGeoTube *tub = new TGeoTube(5 * cm, 6 * cm, 20 * cm);
    TGeoVolume *vtub = new TGeoVolume("Tub", tub, mpb);
    TGeoTube *tubIn = new TGeoTube(0 * cm, 5 * cm, 20 * cm);
    TGeoVolume *vtubIn = new TGeoVolume("TubIn", tubIn, mcl);
    TGeoConeSeg *con = new TGeoConeSeg(10 * cm, 0, 5 * cm, 0, 40 * cm , 45 * degree, (360 - 45)*degree);
    TGeoVolume *vcon = new TGeoVolume("Con", con, mal);
    TGeoBBox *box = new TGeoBBox(10 * cm, 10 * cm, 10 * cm);
    TGeoVolume *vbox = new TGeoVolume("Box", box, mal);
    TGeoPgon *pgon = new TGeoPgon(0, 360, 3, 2);
    {
        pgon->DefineSection(0, 0, 0, 10 * cm);
        pgon->DefineSection(1, 40 * cm, 0, 25 * cm);
    }
    TGeoVolume *vpgon = new TGeoVolume("PGon", pgon, mal);
    TGeoTrap* trap = new TGeoTrap(10*cm,15,30, 10*cm,10*cm,5*cm,0, 10*cm,10*cm,5*cm,0);
    TGeoVolume *vtrap = new TGeoVolume("Trap", trap, mal);

    TGeoHMatrix *mpcon1 = new TGeoHMatrix(); mpcon1->SetDz(20*cm);
    TGeoHMatrix *mpcon2 = new TGeoHMatrix(); mpcon2->SetDz(70*cm);
    vtub->AddNode(vtubIn,0);
    vpcon->AddNode(vtub,0,mpcon1);
    vpcon->AddNode(vsph,1,mpcon2);
    TGeoHMatrix *mtop1 = new TGeoHMatrix(); mtop1->RotateY(-90); mtop1->SetDx(100* cm);
    TGeoHMatrix *mtop2 = new TGeoHMatrix(); mtop2->RotateZ(180); mtop2->RotateY(90); mtop2->SetDx(-10 * cm);
    TGeoHMatrix *mtop3 = new TGeoHMatrix(); mtop3->SetDx(-30* cm);
    TGeoHMatrix *mtop4 = new TGeoHMatrix(); mtop4->SetDx(-50* cm);
    TGeoHMatrix *mtop5 = new TGeoHMatrix(); mtop5->RotateY(-90); mtop5->SetDx(-60* cm);
    TGeoHMatrix *mtop6 = new TGeoHMatrix(); mtop6->RotateY(-90); mtop6->SetDx(-40* cm); mtop6->SetDy(30* cm);
    top->AddNode(vpcon,0,mtop1);
    top->AddNode(vcon,1,mtop2);
    top->AddNode(vbox,2,mtop3);
    top->AddNode(vbox,3,mtop4);
    top->AddNode(vpgon,4,mtop5);
    top->AddNode(vtrap,5, mtop6);
    

    gGeoManager->SetTopVolume(top);
    gGeoManager->CloseGeometry();
    top->SetLineColor(kMagenta);
    gGeoManager->SetTopVisible();

    gGeoManager->GetTopVolume()->Draw();

    gSystem->Load("libTShield.so");
    TShield b;
    b.SetGeometry(gGeoManager);
    b.PrintGeometry();
}
