#include "Rays.h"
#include "FileHelper.h"
#include <stdexcept>
#include <iostream>

using std::exception;
using std::string;
using std::runtime_error;

namespace TpcAlignmentLaserRays {
Rays::Rays(Debug const &inDebug) : vDebug(inDebug) {}

void Rays::AddTrack(Track const &t)
{
   vRays.push_back(t);
}
void Rays::AddTrack(Track &&t)
{
   vRays.push_back(t);
}
RaysType Rays::Get()
{
   return vRays;
}
RaysType const &Rays::Get() const
{
   return vRays;
}
Track const &Rays::GetTrack(size_t i) const
{
   return vRays[i];
}

void Rays::ReadRays(const char *RaysFilename, FILE *fp)
{
   int    vTrackNumber{0};
   char   vEOL;
   string errorMsg(RaysFilename);
   errorMsg.append(": FillInRays: ");
   int n{0};
   while (!feof(fp)) {
      ++vTrackNumber;
      Track vTrack(vDebug, vTrackNumber);
      vTrack.ReadTrack(fp, errorMsg);
      vRays.push_back(vTrack);
      n = fscanf(fp, "%c", &vEOL);
      if (vEOL != '\n') {
         vDebug.Print("Wrong EOL char %c: ", vEOL);
      }
      vDebug.Print("track close\n");
   }
}

FILE *Rays::GetFileStream(const char *RaysFilename, string mode)
{
   FILE *fp;
   try {
      fp = FileHelper::get_file_descriptor(RaysFilename, "r");
   } catch (exception const &ex) {
      throw runtime_error(
         string("Read4File: ").append(RaysFilename).append(" : can't open file. Error: ").append(ex.what()).c_str());
   } catch (...) {
      throw runtime_error(string(RaysFilename).append(": can't open file ").c_str());
   }
   return fp;
}

void Rays::Read4File(const char *RaysFilename)
{
   string mode{"r"};
   FILE  *fp = GetFileStream(RaysFilename, mode);
   try {
      ReadRays(RaysFilename, fp);
   } catch (exception const &ex) {
      fclose(fp);
      throw runtime_error(
         string("Error while reading from the file: ").append(RaysFilename).append(": ").append(ex.what()).c_str());
   } catch (...) {
      fclose(fp);
      throw runtime_error(string("Unknown error while reading from file: ").append(RaysFilename).c_str());
   }
   fclose(fp);
}

void Rays::Write2File(const char *RaysFilename)
{
   FILE *fp;
   try {
      fp = FileHelper::get_file_descriptor(RaysFilename, "w");
   } catch (exception const &ex) {
      throw runtime_error(string(RaysFilename).append(": can't open file. Error: ").append(ex.what()));

   } catch (...) {
      throw runtime_error(string(RaysFilename).append(": can't open file ").c_str());
   }
   string errMsg;
   try {
      SaveRays(fp, RaysFilename);
   } catch (exception const &ex) {
      errMsg.append("Error while saving Rays to the file: ").append(RaysFilename).append(": ").append(ex.what());
   } catch (...) {
      errMsg.append("Unknown error while saving Rays to the file: ").append(RaysFilename);
   }

   fclose(fp);
   if (!errMsg.empty()) {
      throw runtime_error(errMsg);
   }
}

void Rays::SaveRays(FILE *fp, const char *filename)
{
   for (Track const &track : vRays) {
      track.WriteTrack(fp, string("Save Rays to the file: ").append(filename));
   }
}
} // namespace TpcAlignmentLaserRays
