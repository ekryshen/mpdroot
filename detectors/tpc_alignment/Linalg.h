#ifndef LINALG_HH
#define LINALG_HH

namespace TpcAlignment{

////////////////////  Алгебраические программы  ////////////////////
//        Матрицы размещаются по столбцам
class Linalg
{
public:
Linalg() = delete;
/**
 * @brief Обнуление вектора
 * 
 * @param A - вектор
 * @param P - размерность
 */
static void clear(double* A,int P) ;

/**
 * @brief Копирование вектора
 * 
 * @param A 
 * @param B 
 * @param P 
 */
static void copy(double* A,double* B,int P);

/**
 * @brief Линейная комбинация векторов C=a*A+b*B
 * 
 * @param a 
 * @param A 
 * @param b 
 * @param B 
 * @param n 
 * @param C 
 */
static void lcomb(double a,double *A,double b,double *B,int n,double *C);
/**
 * @brief Транспонирование матрицы 
 * 
 * @param A 
 * @param B 
 * @param P 
 * @param Q 
 */
static void transpose(double* A,double* B,int P,int Q);

/**
 * @brief Умножение матрицы A на матрицу B, результат C
 * 
 * @param A 
 * @param B 
 * @param C 
 * @param P 
 * @param Q 
 * @param R 
 */
static void mult(double* A,double* B,double* C,int P,int Q,int R);

/**
 * @brief  Обращение матрицы A методом Гаусса
 * 
 * @param A - матрица
 * @param Ainv - инвариант
 * @param P - размерность
 */
static void inverse(double* A,double* Ainv,int P);

/**
 * @brief Скалярное произведение векторов
 * 
 * @param A  - вектор
 * @param B - вектор
 * @param P - размерность
 * @return double 
 */
static double cprod(double* A,double* B,int P);

/**
 * @brief Косинус излома
 * 
 * @param p1 
 * @param p2 
 * @param p3 
 * @param n 
 * @return double 
 */
static double cosangle(double* p1,double* p2,double* p3, int n);

/**
 * @brief Векторное произведение A*B = C
 * 
 * @param A - вектор
 * @param B - вектор
 * @param C - Результирующий вектор
 */
static void vprod(double* A,double* B,double *C);

/**
 * @brief Нормировка вектора
 * 
 * @param A 
 * @param P 
 */
static void Normalize(double *A,int P) ;

};

}

#endif