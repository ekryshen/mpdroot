TChain* Chain(Int_t nFiles, TString firstFile)
{
  // Get first file number
  Int_t leng = firstFile.Length(), i1 = 0, i2 = 0;
  //cout << leng << endl;
  TString numb, prefix, suffix, symb, symb0;
  //cout << numb.Length() << endl;
  for (Int_t i = leng-1; i > -1; --i) {
    symb = TString(firstFile(i,1));
    if (symb == "_") {
      prefix = firstFile(0,i+1);
      i1 = i + 1;
      break;
    } else if (symb == ".") {
      suffix = firstFile(i,leng-i);
      i2 = i - 1;
    }
  }
  numb = TString(firstFile(i1,i2-i1+1));

  Int_t numb0 = numb.Atoi();
  cout << numb << endl;
  cout << numb0 << endl;
  cout << prefix << endl;
  cout << suffix << endl;

  TChain *chain = new TChain("mpdsim");
  TString fileName;
  nFiles += numb0;
  for (Int_t i = numb0; i < nFiles; ++i) {
    fileName = prefix;
    fileName += i;
    fileName += suffix;
    if (!gSystem->FindFile("./",fileName)) break;
    chain->AddFile(fileName);
  }
  chain->ls();
  return chain;
}
