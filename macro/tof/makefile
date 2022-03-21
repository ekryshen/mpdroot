read_MC_DST: 
	g++ -g3 -o $@ read_MC_DST.cxx  -ffriend-injection -I`root-config --incdir` -I$(PLUTO_INCLUDE_DIR) -L$(PLUTO_LIBRARY_DIR) -I$(VMCWORKDIR)/simulation/mcStack -I$(VMCWORKDIR)/field -I$(VMCWORKDIR)/simulation/generators -I$(VMCWORKDIR)/core/mpdPassive -I$(VMCWORKDIR)/trackbase -I$(VMCWORKDIR)/detectors/bbc -I$(VMCWORKDIR)/detectors/ffd -I$(VMCWORKDIR)/reconstruction/tracking/kalman -I$(VMCWORKDIR)/reconstruction/tracking/lheTrack -I$(VMCWORKDIR)/detectors/sts -I$(VMCWORKDIR)/detectors/etof -I$(VMCWORKDIR)/detectors/tpc -I$(VMCWORKDIR)/detectors/tof -I$(VMCWORKDIR)/detectors/zdc -L`root-config --libdir` -L$(VMCWORKDIR)/build/lib -lBbc -lCpc -lEmc -lEtof -lFfd -lFsa -lSft -lStrawendcap -lSts -lStt -ltpc -lZdc -lTof -lField -lPassive -lKalman -lLHETrack -lMpdData -lMpdgenerators -lPluto -lMpdBase -lBase -lMCStack -lTrkBase -lParBase -lGeoBase -lGeane -lGen -lXMLIO -lMinuit -lGeom -lGeomPainter -lVMC -lEG `root-config --glibs`


clean:
	@rm -f read_MC_DST
