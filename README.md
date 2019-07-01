# Photometry-tools

A python script to do photometry on astronomical Fits files

The process is:
	Reduction; using calibration files provided and writing calibrated files to disk
	
	Astrometry: runs plate solving using Dan Pereley's AutoAstrometry script and updates WCS headers to the files

	Photometry: Uses Sextractor and PSFex to identify sources in each image and writes the table to a file.
