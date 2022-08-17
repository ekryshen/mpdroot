#include "Linalg.h"
#include <cmath>

using namespace std;

namespace TpcAlignment {

	void Linalg::clear(double* A, int P) {
		for (int i = 0; i < P; i++) A[i] = 0.;
	}

	void Linalg::copy(double* A, double* B, int P) {
		for (int i = 0; i < P; i++) B[i] = A[i];
	}

	void Linalg::lcomb(double a, double* A, double b, double* B, int n, double* C) {
		for (int i = 0; i < n; i++) C[i] = a * A[i] + b * B[i];
	}

	void Linalg::transpose(double* A, double* B, int P, int Q) {
		// printf("B = %d, P = %d\n",P,Q);
		for (int i = 0; i < P; i++) {
			for (int j = 0; j < Q; j++) {
				B[j + i * Q] = A[i + j * P];
			}
		}
	}

	void Linalg::mult(double* A, double* B, double* C, int P, int Q, int R) {
		double r;                 					 //printf("B = %d, P = %d, R = %d\n",P,Q,R);
		for (int i = 0; i < P; i++)
			for (int j = 0; j < R; j++) {
				r = 0;
				for (int k = 0; k < Q; k++)  r = r + A[i + k * P] * B[k + j * Q];             //printf("i = %d   j = %d   k = %d\n",i,j,k); }
				C[i + j * P] = r;
			}
	}

	void Linalg::inverse(double* A, double* Ainv, int P) {
		int i, j, k; double r;
		for (i = 0; i < P * P; i++) Ainv[i] = 0;
		for (i = 0; i < P; i++) Ainv[i + i * P] = 1;
		for (i = 0; i < P; i++) {
			r = A[i + i * P];
			for (j = i; j < P; j++) A[i + j * P] = A[i + j * P] / r;
			for (j = 0; j < P; j++) Ainv[i + j * P] = Ainv[i + j * P] / r;
			for (k = i + 1; k < P; k++) {
				r = A[k + i * P];
				for (j = i; j < P; j++) A[k + j * P] = A[k + j * P] - r * A[i + j * P];
				for (j = 0; j < P; j++) Ainv[k + j * P] = Ainv[k + j * P] - r * Ainv[i + j * P];
			}
		}
		for (i = P - 2; i >= 0; i--) {
			for (j = 0; j < P; j++)
				for (k = i + 1; k < P; k++)
					Ainv[i + j * P] = Ainv[i + j * P] - A[i + k * P] * Ainv[k + j * P];
		}
	}

	double Linalg::cprod(double* A, double* B, int P) {
		double r = 0;
		for (int i = 0; i < P; i++) {
			r = r + B[i] * A[i];
		}
		return r;
	}

	double Linalg::cosangle(double* p1, double* p2, double* p3, int n) {
		double v1[3], v2[3];
		lcomb(1., p2, -1., p1, n, v1);
		lcomb(1., p3, -1., p2, n, v2);
		return cprod(v1, v2, n) / sqrt(cprod(v1, v1, n) * cprod(v2, v2, n));
	}


	void Linalg::vprod(double* A, double* B, double* C) {
		C[0] = A[1] * B[2] - A[2] * B[1];
		C[1] = A[2] * B[0] - A[0] * B[2];
		C[2] = A[0] * B[1] - A[1] * B[0];
	}
	//FIXME: Never used
	void Linalg::Normalize(double* A, int P) {
		double r = sqrt(cprod(A, A, P));
		//FIXME: sqrt < 0 ? fabs(r)<tolerance?(1.E-32)?
		if (r > 0)
		{
			for (int i = 0; i < P; i++) {
				A[i] = A[i] / r;
			}
		}
	}
}