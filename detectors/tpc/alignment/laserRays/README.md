<H1>
Library entry point is a Runner class.
</H1>

<b>
There is two mode to run the program. <br>
The first one is calibration:<br>
</b>



	// Init Runner:
	Runner R;

	// Set Model data
	R.SetModelData("LaserRays.txt");

	// Set default Correction Coefficient Matrix
	R.SetCorrectionMatrix({}, false);

	// Calibrate
	const int vNumberOfCalibrationIterations{ 6 };
	R.Calibrate(vNumberOfCalibrationIterations);

	// Save Correction Coeff matrix 
	R.SaveAMR2Files("A.out", "MR.out");




<br> 
<b>
Then you may correct input traks: 
</b>

    R.SetCorrectionMatrix("A.out", true); 
    R.CorrectTracks(inputTracksFile, outputTracsFile);  


