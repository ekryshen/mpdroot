#ifndef COMMON_FUNCTIONS_SOURCE_INCLUDED
#define COMMON_FUNCTIONS_SOURCE_INCLUDED

// check whether file exists
bool CheckFileExist(TString &fileName)
{
   gSystem->ExpandPathName(fileName);
   if (gSystem->AccessPathName(fileName.Data()) == true) {
      std::cerr << "\nMissing file: " << fileName << "\n";
      return false;
   }
   return true;
}
#endif // #ifndef COMMON_FUNCTIONS_SOURCE_INCLUDED
