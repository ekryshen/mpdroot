#include "FileHelper.h"
#include <iostream>
#include <string>
#include <exception>
#include <tuple>
#include "Track.h"
#include "Enums.h"
#include <fstream>
#include <iomanip>

using std::string;
using std::runtime_error;
using std::to_string;

namespace TpcAlignmentLaserRays {
const string testPrefix{"../../TpcAlignmentMTest/"};
const string releasePrefix{"../../TpcAlignmentRun/"};
const string OutputPath{"output/"};
const string InputPath{"input/"};
const string testSuffix{"_4Test.txt"};
const string releaseSuffix{""};
const string extension{".txt"};

void FileHelper::FileExists(FILE *file, string ErrMsg)
{
   if (file == nullptr) { // NULL
      throw runtime_error(ErrMsg);
   }
}

std::tuple<string, string, string> FileHelper::Path(Solution solution, Direction direction)
{
   string path   = (direction == Direction::input) ? InputPath : OutputPath;
   string prefix = (solution == Solution::test) ? testPrefix : releasePrefix;
   string suffix = (solution == Solution::test) ? testSuffix : releaseSuffix;
   return make_tuple(prefix, path, suffix);
}

string FileHelper::BuildFilePath(Solution solution, Direction direction, string filename)
{
   string prefix, path, suffix;
   std::tie(prefix, path, suffix) = Path(solution, direction);
   return string(prefix).append(path).append("/").append(filename).append(suffix);
}

string FileHelper::BuildFilePath(Solution solution, Direction direction, string filename, string innerPath)
{
   string prefix, path, suffix;
   std::tie(prefix, path, suffix) = Path(solution, direction);
   string inner                   = (innerPath.empty() ? "" : innerPath.append("/"));
   return string(prefix).append(path).append(inner).append(filename).append(suffix);
}

string FileHelper::BuildFilePath(Solution solution, Direction direction, unsigned track, string filename,
                                 string innerPath)
{
   string prefix, path, suffix;
   std::tie(prefix, path, suffix) = Path(solution, direction);
   string filepath{string(prefix).append(path)};

   if (!innerPath.empty()) {
      filepath.append(innerPath).append(to_string(track)).append("/");
   }
   filepath.append(filename).append("_").append(to_string(track)).append(suffix);
   return filepath;
}

string FileHelper::BuildFilePath(Solution solution, Direction direction, string filename, string innerPath,
                                 unsigned track)
{
   string prefix, path, suffix;
   std::tie(prefix, path, suffix) = Path(solution, direction);
   string filepath{string(prefix).append(path)};
   if (!innerPath.empty()) {
      filepath.append(innerPath).append("/");
   }
   filepath.append(filename).append("_").append(to_string(track)).append(suffix);
   return filepath;
}

string FileHelper::BuildFilePath(Solution solution, Direction direction, unsigned track, size_t point, string filename,
                                 string innerPath)
{
   string prefix, path, suffix;
   std::tie(prefix, path, suffix) = Path(solution, direction);
   string filepath{string(prefix).append(path)};
   if (!innerPath.empty()) {
      filepath.append(innerPath).append(to_string(track)).append("/");
   }
   filepath.append(filename).append("_").append(to_string(track)).append("_").append(to_string(point));
   filepath.append(suffix);
   return filepath;
}

FILE *FileHelper::get_file_descriptor(const char *Filename, string mode)
{
   FILE *fp;
   fp = fopen(Filename, mode.c_str());
   try {
      FileHelper::FileExists(fp, string("get_file_descriptor: can't open file").append(Filename));
   } catch (std::exception const &ex) {
      throw ex;
   } catch (...) {
      throw runtime_error(string("Unknown error while opening file: ").append(Filename));
   }
   return fp;
}

void FileHelper::ReadMR_obsolete(string filename, SectorMatrix &vM, SectorMatrix &vR)
{
   FILE *fp = fopen(filename.c_str(), "r");
   if (fp == NULL) {
      throw runtime_error(string("Couldn't open file: ").append(filename));
   }
   for (int s = 0; s < SECTORS; ++s) {
      vM[s].resize(LL);
      vR[s].resize(L_);
      for (size_t i = 0; i < L_; i++) {
         for (size_t j = 0; j < L_; j++) {
            fscanf(fp, "%lf", &vM[s][i + L_ * j]);
         }
         fscanf(fp, "%lf\n", &vR[s][i]);
      }
      fscanf(fp, "\n");
   }
   fclose(fp);
}

void FileHelper::SaveMR(string filename, SectorMatrix const &vM, SectorMatrix const &vR, int precision)
{
   std::ofstream fs(filename.c_str()); // , std::ofstream::binary
   if (!fs.is_open()) {
      throw runtime_error(string("Can't open a file: ").append(filename).append(" to save a MR-coefficients"));
   }

   for (int s = 0; s < SECTORS; s++) {
      for (size_t i = 0; i < L_; i++) {
         for (size_t j = 0; j < L_; j++) {
            fs << std::fixed << std::setprecision(precision) << vM.at(s).at(i + L_ * j) << " ";
         }
         fs << std::fixed << std::setprecision(precision) << vR.at(s).at(i) << " ";
      }
      fs << '\n';
   }
   fs.close();
}

void FileHelper::LoadMR(string filename, SectorMatrix &vM, SectorMatrix &vR)
{
   std::ifstream fs(filename.c_str()); // , std::ofstream::binary
   if (!fs.is_open()) {
      throw runtime_error(string("Can't open a file: ").append(filename).append(" to save a CorrectionMatrix"));
   }

   for (int s = 0; s < SECTORS; ++s) {
      vM[s].resize(LL);
      vR[s].resize(L_);
      for (size_t i = 0; i < L_; i++) {
         for (size_t j = 0; j < L_; j++) {
            fs >> vM[s][i + L_ * j];
         }
         fs >> vR[s][i];
      }
   }
   fs.close();
}

void FileHelper::SaveCorrectionMatrix(std::string filename, SectorMatrix const &vCorrectionMatrix, int precision)
{
   std::ofstream fs(filename.c_str()); // , std::ofstream::binary
   if (!fs.is_open()) {
      throw runtime_error(std::string("Can't open a file: ").append(filename).append(" to savel a CorrectionMatrix"));
   }

   for (size_t k = 0; k < SECTORS; k++) {
      for (size_t i = 0; i < 3; i++) {
         for (size_t j = 0; j < M_; j++) {
            fs << std::fixed << std::setprecision(precision) << vCorrectionMatrix[k][i + 3 * j] << " ";
         }
         fs << '\n';
      }
      fs << '\n';
   }
   fs.close();
}

void FileHelper::LoadCorrectionMatrix(std::string filename, SectorMatrix &vCorrectionMatrix)
{
   std::ifstream fs(filename.c_str()); // , std::ofstream::binary
   if (!fs.is_open()) {
      throw runtime_error(std::string("Can't open a file: ").append(filename).append(" to save a CorrectionMatrix"));
   }

   for (size_t k = 0; k < SECTORS; k++) {
      vCorrectionMatrix.at(k).resize(L_);
      for (size_t i = 0; i < 3; i++) {
         for (size_t j = 0; j < M_; j++) {
            fs >> vCorrectionMatrix[k][i + 3 * j]; // >> fixed >> setprecision(precision)
         }
         //	  '\n'
      }
      //   '\n'
   }
   fs.close();
}
/*******************************************************************************/
void FileHelper::SaveCorrectionMatrix_obsolete(std::string filename, SectorMatrix const &vCorrectionMatrix)
{
   FILE *fp = fopen(filename.c_str(), "w");

   if (fp == NULL) {
      return;
   }

   for (size_t k = 0; k < SECTORS; k++) {
      for (size_t i = 0; i < 3; i++) {
         for (size_t j = 0; j < M_; j++) {
            fprintf(fp, " %f", vCorrectionMatrix[k][i + 3 * j]);
         }
         fprintf(fp, "%c", '\n');
      }
      fprintf(fp, "%c", '\n');
   }
   fclose(fp);
}

void FileHelper::ReadCorrectionMatrix_obsolete(std::string filename, SectorMatrix &vCorrectionMatrix)
{
   FILE *fp{nullptr};
   float r;

   fp = fopen(filename.c_str(), "r");
   if (fp == NULL) {
      throw runtime_error(string("Can't read the file: ").append(filename).append(". Correction matrix is not set!"));
   }
   int readN{0};
   for (size_t k = 0; k < SECTORS; k++) {
      vCorrectionMatrix.at(k).resize(L_);
      for (size_t i = 0; i < 3u; i++) {
         for (size_t j = 0; j < M_; j++) {
            readN = fscanf(fp, " %f", &r);
            if (readN == 1) {
               try {
                  vCorrectionMatrix[k][i + 3u * j] = r;
               } catch (const std::exception &ex) {
                  fclose(fp);
                  throw runtime_error(string("Error readind file :").append(filename).append(" : ").append(ex.what()));
               }
            } else {
               fclose(fp);
               throw runtime_error(string("Error in file: ")
                                      .append(filename)
                                      .append("Line: ")
                                      .append(to_string(k))
                                      .append(", Column: ")
                                      .append(to_string(i * j)));
            }
         }
      }
   }
   fclose(fp);
}

} // namespace TpcAlignmentLaserRays
