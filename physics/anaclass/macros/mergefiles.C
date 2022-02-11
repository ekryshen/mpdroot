void mergefiles(Int_t nStarFile = 340, Int_t nEndFile= 345){

TFileMerger m;

for(Int_t jk=nStarFile; jk<nEndFile; jk++){
char NameOfDir[999];
sprintf(NameOfDir,"/scratch1/maldonado/Soft/Centrality/ana/outputana%d.root",jk);
if(!(m.AddFile(NameOfDir)))continue;
}

m.OutputFile("outputmerged.root");
m.Merge();
}
