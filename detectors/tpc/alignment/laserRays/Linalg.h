#ifndef LINALG_HH
#define LINALG_HH
namespace TpcAlignment {
////////////////////  Linear algebra programmes  ////////////////////
//        Matrices are column oriented
class Linalg {
public:
   Linalg() = delete;
   /**
    * @brief Zeroing of the vector
    *
    * @param A - vector
    * @param P - dimension
    */
   static void clear(double *A, int P);

   /**
    * @brief Copying vector
    *
    * @param A
    * @param B
    * @param P
    */
   static void copy(double *A, double *B, int P);

   /**
    * @brief Linear combination of vectors C=a*A+b*B
    *
    * @param a
    * @param A
    * @param b
    * @param B
    * @param n
    * @param C
    */
   static void lcomb(double a, double *A, double b, double *B, int n, double *C);
   /**
    * @brief Matrix transposing
    *
    * @param A
    * @param B
    * @param P
    * @param Q
    */
   static void transpose(double *A, double *B, int P, int Q);

   /**
    * @brief Multiplying matrix A by matrix B. Result is in C
    *
    * @param A
    * @param B
    * @param C
    * @param P
    * @param Q
    * @param R
    */
   static void mult(double *A, double *B, double *C, int P, int Q, int R);

   /**
    * @brief  Matrix inversion using Gauss method
    *
    * @param A - matrix
    * @param Ainv - inverse
    * @param P - dimension
    */
   static void inverse(double *A, double *Ainv, int P);

   /**
    * @brief Scalar multiplication of vectors
    *
    * @param A  - vector
    * @param B - vector
    * @param P - dimension
    * @return double
    */
   static double cprod(double *A, double *B, int P);

   /**
    * @brief Cosine of vector angle
    *
    * @param p1
    * @param p2
    * @param p3
    * @param n
    * @return double
    */
   static double cosangle(double *p1, double *p2, double *p3, int n);

   /**
    * @brief Vector multiplication of vectors A*B = C
    *
    * @param A - vector
    * @param B - vector
    * @param C - resulting vector
    */
   static void vprod(double *A, double *B, double *C);

   /**
    * @brief Vector normalization
    *
    * @param A
    * @param P
    */
   static void Normalize(double *A, int P);
};

} // namespace TpcAlignment
#endif
