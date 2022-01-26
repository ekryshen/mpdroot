read_MC_DST: 
	g++ -g3 -o $@ read_MC_DST.cxx  -ffriend-injection -I`root-config --incdir` -I$(PLUTO_INCLUDE_DIR) -L$(PLUTO_LIBRARY_DIR) -I$(VMCWORKDIR)/base -I$(VMCWORKDIR)/mcstack -I$(VMCWORKDIR)/parbase -I$(VMCWORKDIR)/field -I$(VMCWORKDIR)/generators -I$(VMCWORKDIR)/passive -I$(VMCWORKDIR)/trackbase -I$(VMCWORKDIR)/detectors/bbc -I$(VMCWORKDIR)/cpc -I$(VMCWORKDIR)/fsa -I$(VMCWORKDIR)/detectors/ffd -I$(VMCWORKDIR)/kalman -I$(VMCWORKDIR)/lhetrack -I$(VMCWORKDIR)/strawendcap -I$(VMCWORKDIR)/sts -I$(VMCWORKDIR)/detectors/etof -I$(VMCWORKDIR)/tpc -I$(VMCWORKDIR)/tof -I$(VMCWORKDIR)/zdc  -I$(VMCWORKDIR)/mpddata -L`root-config --libdir` -L$(VMCWORKDIR)/build/lib -lBbc -lCpc -lEmc -lEtof -lFfd -lFsa -lSft -lStrawendcap -lSts -lStt -ltpc -lZdc -lTof -lField -lPassive -lKalman -lLHETrack -lMpdData -lMpdgenerators -lPluto -lMpdBase -lBase -lMCStack -lTrkBase -lParBase -lGeoBase -lGeane -lGen -lXMLIO -lMinuit -lGeom -lGeomPainter -lVMC -lEG `root-config --glibs`


clean:
	@rm -f read_MC_DST
