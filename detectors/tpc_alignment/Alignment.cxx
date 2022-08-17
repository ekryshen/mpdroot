#include "Alignment.h"
#include "Linalg.h"
#include <cerrno>
#include <cstdio>
#include <iostream>
#include <math.h>
#include "FileHelper.h"
#include "Debug.h"

using namespace std;

#define MAXn 1000
#define PI 3.14159265

namespace TpcAlignment {

	Alignment::Alignment(Debug const & inMode) :
		vDebug(inMode) 
	{
		vDebug.Print("testA start");
	}

	Alignment::~Alignment() { vDebug.Print("testA finished"); }

	void Alignment::SetA2Default() {
		vDebug.Print("SetA2Default\n");
		for (int k = 0; k < 24; k++) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < M_; j++) {
					A[k * L_ + i + 3 * j] = 0.;
				}
				A[k * L_ + i + 3 * i] = 1.;
			}
		}
	}

	void Alignment::LoadA() {
		if (!LoadA4File()) {
			SetA2Default();
			SaveA2File();
		}
	}

	bool Alignment::LoadA4File(const char *filename) {
		FILE *fp{ nullptr };
		float r;
		vDebug.Print("Load A from a file: %s\n", filename);

		fp = fopen(filename, "r");
		if (!FileHelper::FileExists(fp, filename, "can't read file")) {
			return false;
		}
		int readN{ 0 };
		for (int k = 0; k < 24; k++) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < M_; j++) {
					readN = fscanf(fp, " %f", &r);
					if (readN == 1) {
						A[k * L_ + i + 3 * j] = r;
					}
					else {
						cout << "Error in file: " << filename << '\n';
						cout << "Line: " << k << ", Column: " << i * j << '\n';
						fclose(fp);
						return false;
					}
				}
			}
		}
		fclose(fp);
		return true;
	}

	void Alignment::SaveA2File(const char *filename) {
		vDebug.Print("SaveA2File");
		FILE *fp = fopen(filename, "w");

		if (!FileHelper::FileExists(fp, filename, "can't create file"))
		{
			return;
		}

		for (int k = 0; k < 24; k++) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < M_; j++) {
					fprintf(fp, " %f", A[k * L_ + i + 3 * j]);
				}
				fprintf(fp, "%c", '\n');
			}
			fprintf(fp, "%c", '\n');
		}
		fclose(fp);
	}

	void Alignment::ResetCoefficients() {
		vDebug.Print("ResetCoefficients");

		for (int i = 0; i < 24 * LL; i++) {
			M[i] = 0.;
		}
		for (int i = 0; i < 24 * L_; i++) {
			R[i] = 0.;
		}
	}

	void Alignment::Finish() {

		vDebug.Print("Finish");

		//
		//              Решение системы уравнений
		//
		for (int k = 0; k < 24; k++) {
			lsolve(&M[k * LL], &R[k * L_], L_, &A[k * L_]);
		}

		SaveAMRCoeff2File();
		SaveA2File();
	}

	void Alignment::SaveAMRCoeff2File(const char *filename) {
		FILE *fp;
		fp = fopen(filename, "w");
		if (!FileHelper::FileExists(fp, filename, "can't create file"))
		{
			return;
		}

		for (int k = 0; k < 24; k++) {
			for (int i = 0; i < L_; i++) {
				for (int j = 0; j < L_; j++) {
					fprintf(fp, " %f", M[k * LL + i + L_ * j]);
				}
				fprintf(fp, "  %f\n", R[k * L_ + i]);
			}
			fprintf(fp, "%c", '\n');
		}
		fclose(fp);
	}

	void Alignment::SetF(double *f, double x, double y, double z) {
		f[0] = x;
		f[1] = y;
		f[2] = z;
		f[3] = 1;
	}

	void Alignment::AddTrack(int n, double *v) {
		double v2[3] = { 0., 0., 0. };
		double P[3] = { 0., 0., 0. };
		double Q[3] = { 0., 0., 0. };
		double X[3] = { 0., 0., 0. };
		double G[LL];
		double H[L_];
		double r{ 0 };
		double w[3] = { 0., 0., 0. };
		double vv[M_];
		// X=P*t+Q - параметрическое задание трека
		double x, y, z;
		double f[M_], J[9], J1[9];
		int kk[MAXn], k;

		for (int m = 0; m < n; m++) {
			kk[m] = Correct(&v[m * 3 + 0], v2); // k[j] - номер сектора точки
			for (int i = 0; i < 3; i++) {
				Q[i] = Q[i] + v2[i] / n;
			}
			if (m == 0) {
				for (int i = 0; i < 3; i++) {
					P[i] = v2[i];
				}
			}
			if (m == (n - 1)) {
				for (int i = 0; i < 3; i++) {
					P[i] = v2[i] - P[i];
					r = r + P[i] * P[i];
				}
			}
			//	for(int i=0;i<3;i++) v[i+m*3]=v2[i];  				//
			//??????
		}
		r = sqrt(r);
		for (int i = 0; i < 3; i++) {
			P[i] = P[i] / r;
		}
		Linalg::clear(J, 9);
		for (int k = 0; k < n; k++) {
			Linalg::lcomb(1., &v[k * 3], -1., Q, 3, v2);
			for (int i = 0; i < 3; i++) {
				J[i + i * 3] = J[i + i * 3] + Linalg::cprod(v2, v2, 3) - v2[i] * v2[i];
				for (int j = 0; j < 3; j++) {
					if (i == j) {
						J[i + j * 3] = J[i + j * 3] + v2[i] * v2[j];
					}
				}
			}
		}
		Linalg::inverse(J, J1, 3);
		EU3(J1, P);

		vDebug.Print("AddTrack   %f %f %f   %f %f %f\n", P[0], P[1], P[2], Q[0], Q[1],
			Q[2]);

		//
		//  Накопление матрицы системы уравнений
		//
		double angle{ 30. / 180. * PI };
		for (int m = 0; m < n; m++) {
			k = kk[m];
			//FIXME: Next 6 lines?
			w[0] = 80. * sin(k * angle);
			w[1] = 80. * cos(k * angle);
			w[2] = 0;
			for (int i = 0; i < 3; i++) {
				w[i] = 0.;
			}

			vDebug.Print(" Segment %d\n", k);

			for (int i = 0; i < 3; i++) {
				X[i] = v[m * 3 + i];
			}
			for (int i = 0; i < 3; i++) {
				vv[i] = X[i] - w[i];
			}
			vv[3] = 1;
			x = vv[0];
			y = vv[1];
			z = vv[2];
			SetF(f, x, y, z);

#include "G.h"
#include "H.h"
			for (int i = 0; i < L_; i++) {
				for (int j = 0; j < L_; j++) {
					M[k * LL + i + j * L_] = M[k * LL + i + j * L_] + G[i + j * L_];
				}
				R[k * L_ + i] = R[k * L_ + i] - H[i];
			}
		}
	}

	int Alignment::Correct(double *v1, double *v2) {
		double fi;
		double r;
		double w[3] = { 0., 0., 0. };
		double vv[M_];
		double x, y, z;
		double f[M_];
		int k;
		fi = atan2(v1[0], v1[1]);
		k = CalculateK(fi, v1[2]);
		//FIXME: Next 6 lines?
		w[0] = 80. * sin(k * 30. / 180. * PI);
		w[1] = 80. * cos(k * 30. / 180. * PI);
		w[2] = 0;
		for (int i = 0; i < 3; i++) {
			w[i] = 0.;
		}
		//FIXME: Next 3 lines?:  Vi-0=Vi ?
		for (int i = 0; i < 3; i++) {
			vv[i] = v1[i] - w[i];
		}
		vv[3] = 1.;
		x = vv[0];
		y = vv[1];
		z = vv[2];
		SetF(f, x, y, z);
		for (int i = 0; i < 3; i++) {
			r = w[i];
			for (int j = 0; j < M_; j++) {
				r = r + A[k * L_ + i + 3 * j] * vv[j];
			}
			//  ????????
			v2[i] = r;
		}
		return k;
	}

	int Alignment::CalculateK(double fi, double z)
	{
		int k = static_cast<int>(floor(fi / PI * 6 + 12.5));
		k = k % 12;
		if (z < 0.) {
			k = k + 12;
		}
		return k;
	}

	int Alignment::lsolve(double *A, double *B, int n, double *X) {
		double r;
		int m;
		for (int i = 0; i < n; i++) {
			m = 0;
			r = 0.;
			for (int j = 0; j < n; j++) {
				if (abs(A[i + j * n]) > r) {
					m = j;
					r = abs(A[i + j * n]);
				}
			}
			r = A[i + m * n];
			// FIXME: double == 0.?
			if (r == 0.) {
				return 0;
			}
			for (int j = 0; j < n; j++) {
				A[i + j * n] = A[i + j * n] / r;
			}
			B[i] = B[i] / r;
			A[i + m * n] = 1.;
			for (int k = i + 1; k < n; k++) {
				r = A[k + m * n];
				for (int j = 0; j < n; j++) {
					A[k + j * n] = A[k + j * n] - r * A[i + j * n];
				}
				B[k] = B[k] - r * B[i];
				A[k + m * n] = 0.;
			}
		}  // 

		for (int i = n - 1; i >= 0; i--) {
			m = 0;
			for (int j = 0; j < n; j++) {
				// FIXME: == double
				if (A[i + j * n] == 1.) {
					m = j;
				}
			}
			X[m] = B[i];
			for (int k = i - 1; k >= 0; k--) {
				r = A[k + m * n];
				B[k] = B[k] - r * B[i];
				A[k + m * n] = 0.;
			}
		}

		return 1;
	}

	double Alignment::EU3(double *A, double *P) {
		double lambda = 0, lambda1, P1[3];
		do {
			lambda1 = lambda;
			Linalg::mult(A, P, P1, 3, 3, 1);
			lambda = sqrt(Linalg::cprod(P1, P1, 3));
			for (int i = 0; i < 3; i++) {
				P[i] = P1[i] / lambda;
			}
		} while (abs(lambda1 / lambda - 1) > 1.e-10);
		return lambda;
	}

} // namespace TpcAlignment