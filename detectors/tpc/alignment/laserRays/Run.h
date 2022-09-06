#ifndef RUN_HH
#define RUN_HH
namespace TpcAlignment {
class Debug;
class Alignment;

class Run {
private:
   // const char *vTrackFilename;
   // const char *vCoefficientFilename;
   Alignment   &vAlignment;
   Debug const &vDebug;

public:
   Run(Alignment &inAlignment, Debug const &inDebug);
   ~Run();
   /**
    * @brief Калибровка: подача измеренных точек track из файла
    *
    * @param TrackFilename
    */
   void Calibration(const char *TrackFilename);
   /**
    * @brief подача измеренных точек из файла и запись скорректированных точек в
    * файл
    *
    * @param f1 - inputfile
    * @param f2 - outputfile
    */
   void Measurement(const char *f1, const char *f2);
};
} // namespace TpcAlignment
#endif
