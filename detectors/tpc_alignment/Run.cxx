#include "Run.h"
#include "FileHelper.h"
#include <cstdio>
#include <exception>
#include <iostream>
#include <string>
#include "Alignment.h"
#include "Debug.h"

#define MAXn 1000

using namespace std;

namespace TpcAlignment {
	Run::Run(Alignment& inAlignment, Debug const& inDebug) :
		vAlignment(inAlignment),
		vDebug(inDebug)
	{}
	Run::~Run() {}

	void Run::Calibration(const char* TrackFilename) {

		vDebug.Print("Calibration started\n");
		double v[MAXn * 3];

		vAlignment.ResetCoefficients();

		FILE* fp;
		float r;
		fp = fopen(TrackFilename, "r");
		if (!FileHelper::FileExists(fp, TrackFilename, "can't open file")) {
			return;
		}
		int PointsOnTrack{ 0 };
		int readNumber{ 0 };
		int vTrackNumber{ 0 };
		char vEOL;

		while (!feof(fp)) {
			PointsOnTrack =
				FileHelper::ReadPointsOnTrack(fp, TrackFilename, "Calibration");
			if (PointsOnTrack == 0) {
				break;
			}
			vDebug.Print("track open\n");
			++vTrackNumber;
			for (int k = 0; k < PointsOnTrack; k++) {
				for (int i = 0; i < 3; i++) {
					readNumber = fscanf(fp, "%f", &r);
					if (readNumber == 0) {
						const char* coordinate =
							(i == 0 ? "X"
								: (i == 1 ? "Y" : (i == 2 ? "Z" : "Unknown coordinate")));
						cout << TrackFilename << " : "
							<< "Wrong number on line: " << vTrackNumber << ", point: " << k
							<< ", coordinate: " << coordinate << '\n';
						fclose(fp);
						return;
					}
					v[k * 3 + i] = r;
				}
				vDebug.Print("get point %f %f %f\n", v[k * 3 + 0], v[k * 3 + 1],
					v[k * 3 + 2]);
			}
			vAlignment.AddTrack(PointsOnTrack, v);
			readNumber = fscanf(fp, "%c", &vEOL);
			if(vEOL != '\n')
			{
				vDebug.Print("Wrong EOL char %c: ", vEOL);
			}
			vDebug.Print("track close\n");
		}
		fclose(fp);
		vAlignment.Finish();
		vDebug.Print("Calibration finished\n");
	}

	void Run::Measurement(const char* inputName, const char* outputName) {
		vDebug.Print("Measurement started\n");
		FILE* inputStream, * outputStream;
		int TrackNumber{ 0 };
		int PointsOnTrack{ 0 };
		double v1[3] = { 1., 0., 0. };
		double v2[3]; //  ????????
		char eol;

		inputStream = fopen(inputName, "r");
		if (!FileHelper::FileExists(inputStream, inputName, "can't open file")) {
			return;
		}
		outputStream = fopen(outputName, "w");
		if (!FileHelper::FileExists(outputStream, outputName, "can't open file")) {
			return;
		}

		while (!feof(inputStream)) {
			PointsOnTrack =
				FileHelper::ReadPointsOnTrack(inputStream, inputName, "Measurement");
			if (PointsOnTrack == 0) {
				break;
			}
			if (!FileHelper::WritePointsOnTrack(outputStream, PointsOnTrack, outputName,
				"Measurement")) {
				break;
			}
			++TrackNumber;
			vDebug.Print("%d track is opened\n", TrackNumber);

			for (int Point = 0; Point < PointsOnTrack; Point++) {
				if (!FileHelper::ReadXYZ(v1[0], v1[1], v1[2], inputStream, inputName,
					"Measurement", TrackNumber, Point)) {
					break;
				}

				vDebug.Print("Point(x,y,z): (%f, %f, %f)\n", v1[0], v1[1], v1[2]);
				vAlignment.Correct(v1, v2);
				if (!FileHelper::WriteXYZ(v2[0], v2[1], v2[2], outputStream, outputName,
					"Measurement", TrackNumber, Point)) {
					break;
				}
			}

			fscanf(inputStream, "%c", &eol);
			if (eol != '\n')
			{
				vDebug.Print("It's not a end of line: %c", eol);
				break;
			}
			fprintf(outputStream, "\n");
			vDebug.Print("%d track is closed\n", TrackNumber);
		}
		fclose(inputStream);
		fclose(outputStream);
		vDebug.Print("Measurement finished\n");
	}
} // namespace TpcAlignment