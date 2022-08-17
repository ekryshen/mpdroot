#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH

#include <cstdio>
#define M_ 4
#define N_ 3
#define L_ 3 * M_
#define LL L_ *L_

namespace TpcAlignment {
	class Debug;
/**
 * @brief Процедуры коррекции измеренных точек
 *
 */
class Alignment {
private:
  double A[24 * L_];
  double M[24 * LL];
  double R[24 * L_];
  const char *vInputFilename;
  const char *vOutputFilename;
  Debug const &vDebug;

public:
  Alignment(Debug const & inMode);
  ~Alignment();
  /**
   * @brief Обнуление коеффициентов системы уравнений
   *
   */
  void ResetCoefficients();

  /**
  * @brief Load A from file or set default values
  *
  */
  void LoadA();

  void Finish();

  /**
  * @brief X=P*t+Q - параметрическое задание трека
  *
  * @param n - размерность
  * @param v - массив
  */
  void AddTrack(int n, double *v);
  int Correct(double *v1, double *v2);
protected:
  /**
  * @brief Load A from a file
  *
  */
  bool LoadA4File(const char *filename = "testA.dat");
  /**
  * @brief Set A to default values
  *
  */
  void SetA2Default();


  void SaveA2File(const char *filename = "testA.dat");
 
  void SaveAMRCoeff2File(const char *filename = "outAMR.dat");

  /**
   * @brief Решение системы уравнений
   *
   * @param A
   * @param B
   * @param n
   * @param X
   * @return int
   */
  static int lsolve(double *A, double *B, int n, double *X);
  static double EU3(double *A, double *P); 
  /**
  *@brief 
  *
  */
  static void SetF(double *f, double x, double y, double z);
  static int CalculateK(double fi, double z);
};

} // namespace TpcAlignment

#endif