#!/usr/bin/env python
#Filename: 2HDECAY.py


 ##########################################################################################################################################
#                                                                                                                                          #
#                                                                2HDECAY                                                                   #
#                                                                                                                                          #
#   Purpose:    A program for the calculation One-Loop Electroweak Corrections to Higgs Decays in the                                      #
#               Two-Higgs-Doublet Model (2HDM) Including State-of-the-Art QCD Corrections                                                  #
#   Authors:    Dr. Marcel Krause (marcel.krause@alumni.kit.edu)                                                                           #
#               Prof. Dr. M. Margarete Muehlleitner (margarete.muehlleitner@kit.edu)                                                       #
#               Dr. Michael Spira (michael.spira@psi.ch)                                                                                   #
#   Version:    1.1.4                                                                                                                      #
#   Date:       04.11.2019                                                                                                                 #
#   Copyright:  Copyright (C) 2018-2019, Marcel Krause, Milada Margarete Muehlleitner and Michael Spira                                    #
#   License:    GNU General Public License (GNU GPL-3.0-or-later)                                                                          #
#                                                                                                                                          #
#               2HDECAY is released under GNU General Public License (GNU GPL-3.0-or-later).                                               #
#               This program is free software: you can redistribute it and/or modify it under the terms of the                             #
#               GNU General Public License as published by the Free Software Foundation, either version 3 of                               #
#               the License, or any later version.                                                                                         #
#                                                                                                                                          #
#               This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;                                  #
#               without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                  #
#               See the GNU General Public License for more details.                                                                       #
#                                                                                                                                          #
#               You have received a copy (LICENSE.md) of the GNU General Public License along with this program                            #
#               in the 2HDECAY root directory.                                                                                             #
#                                                                                                                                          #
#   Citation:   When you use this program, please acknowledge the work of our and other groups by citing the                               #
#               following papers:                                                                                                          #
#                   The manual for 2HDECAY:                                                                                                #
#                    - M. Krause, M. Muhlleitner, M. Spira, Comp. Phys. Commun. 246 (2020) 106852, arXiv:1810.00768 (hep-ph)               #
#                   The manuals for HDECAY:                                                                                                #
#                    - A. Djouadi, J. Kalinowski, M. Spira, Comp. Phys. Commun. 108 (1998) 56, hep-ph/9704448                              #
#                    - A. Djouadi, J. Kalinowski, M. Muhlleitner, M. Spira, arXiv:1801.09506 (hep-ph)                                      #
#                   The papers on the EW correction to the 2HDM decays:                                                                    #
#                    - M. Krause, R. Lorenz, M. Muhlleitner, R. Santos, H. Ziesche, JHEP 1609 (2016) 143, arXiv:1605.04853 (hep-ph)        #
#                    - M. Krause, M. Muhlleitner, R. Santos, H. Ziesche, Phys.Rev. D95 (2017) no.7, 075019, arXiv:1609.04185 (hep-ph)      #
#                    - A. Denner, S. Dittmaier, J.-N. Lang, J. High Energ. Phys. (2018) 2018: 104, arXiv:1808.03466 (hep-ph)               #
#                   The publication of LoopTools:                                                                                          #
#                    - T. Hahn, M. Perez-Victoria, Comp. Phys. Commun. 118 (1999) 153-165, hep-ph/9807565                                  #
#                                                                                                                                          #
 ##########################################################################################################################################


#------------------------------#
#         Import Modules       #
#------------------------------#
import sys
import os
from shutil import copyfile, rmtree
from math import pi, sqrt
import subprocess
# import multiprocessing
import CommonFunctions      # Provides common, often used functions for different scripts of 2HDECAY

#-------------------------#
#        Settings         #
#-------------------------#
# WARNING: do not change these settings if you do not know what they do!
lineToInsert = 124      # This is the line at which the temporary input file ends and at which we append the electroweak corrections
lineWhereAlphaAtMZ = 26 # This is the line at which in the temporary input file the fine-structure constant at the Z boson mass MZ is specified
lineWhereGFCalc = 28    # This is the line at which in the temporary input file the calculated Fermi constant GFCALC is specified
lineWhereMZ = 31        # This is the line at which in the temporary input file the Z boson mass MZ is specified
lineWhereMW = 32        # This is the line at which in the temporary input file the W boson mass MW is specified
lineWhereOSMC = 22      # This is the line at which the OS MC value has to be inserted in the temporary input file
lineWhereOSMB = 23      # This is the line at which the OS MB value has to be inserted in the temporary input file
lineWhereParamType = 56 # This is the line at which the parameter type is specified
lineWhereRefScheme = 59 # This is the line at which the reference renormalization scheme is specified

#----------------------------#
#        Main Program        #
#----------------------------#

if __name__ == "__main__":
	# Print the welcome screen
	print('''
+---------------------------------------+
|                                       |
|             2HDECAY 1.1.4             |
|                                       |
|                             /         |
|                            /          |
|                           /           |
|                      --- /            |
|      ______________/     \            |
|                    \     /            |
|                      --- \            |
|                           \           |
|                            \          |
|                             \         |
|                                       |
+---------------------------------------+

When you use this program please cite:
	The manual for 2HDECAY:
	 - M. Krause, M. Muhlleitner, M. Spira, Comp. Phys. Commun. 246 (2020) 106852, arXiv:1810.00768 (hep-ph)
	The manuals for HDECAY:
	 - A. Djouadi, J. Kalinowski, M. Spira, Comp. Phys. Commun. 108 (1998) 56, hep-ph/9704448
	 - A. Djouadi, J. Kalinowski, M. Muhlleitner, M. Spira, arXiv:1801.09506 (hep-ph)
	The papers on the EW correction to the 2HDM decays:
	 - M. Krause, R. Lorenz, M. Muhlleitner, R. Santos, H. Ziesche, JHEP 1609 (2016) 143, arXiv:1605.04853 (hep-ph)
	 - M. Krause, M. Muhlleitner, R. Santos, H. Ziesche, Phys.Rev. D95 (2017) no.7, 075019, arXiv:1609.04185 (hep-ph)
	 - A. Denner, S. Dittmaier, J.-N. Lang, J. High Energ. Phys. (2018) 2018: 104, arXiv:1808.03466 (hep-ph)
	The publication of LoopTools:
	 - T. Hahn, M. Perez-Victoria, Comp. Phys. Commun. 118 (1999) 153-165, hep-ph/9807565

2HDECAY is released under GNU General Public License (GNU GPL-3.0-or-later). This program is free software: 
you can redistribute it and/or modify it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
for more details.

You have received a copy (LICENSE.md) of the GNU General Public License along with this program in the 2HDECAY
root directory.

Copyright 2018-2019, Marcel Krause, Milada Margarete Muehlleitner and Michael Spira.
	''')

	# Get a list of all input files
	inputPath = "Input"
	inputFileList = os.listdir(inputPath)
	if '..' in inputFileList:
		inputFileList.remove('..')
	if '.' in inputFileList:
		inputFileList.remove('.')

	# Iterate over all input files
	for inputFileTemp in inputFileList:
		print("Calculating corrections for input file " + inputFileTemp + " ...\n")

		# Copy the input file to the HDECAY subfolder
		print("Copying input files into HDECAY folder...")
		filenameIn = "Input" + os.sep + inputFileTemp
		filenameOut = "HDECAY" + os.sep + "hdecay.in"
		# Remove any existing input and fermion masses file in HDECAY and prepare a template mixing angles file
		if os.path.isfile(filenameOut):
			os.remove(filenameOut)
		if os.path.isfile("HDECAY" + os.sep + "fermionmasses.dat"):
			os.remove("HDECAY" + os.sep + "fermionmasses.dat")
		if os.path.isfile("HDECAY" + os.sep + "alphaandbeta.dat"):
			os.remove("HDECAY" + os.sep + "alphaandbeta.dat")
		copyfile(filenameIn, filenameOut)
		with open("HDECAY" + os.sep + "alphaandbeta.dat", 'w') as fileHandler:
			fileHandler.write(' alpha = 0D0\n beta  = 0D0\n')
		print("... done.\n")

		# Let HDECAY run in minimal mode to produce the fermion mass file
		print("Starting HDECAY in minimal mode...")
		os.chdir('HDECAY')
		prompt = ['./run', '1']
		subprocess.call(prompt, stdin=None, stdout=None, stderr=None, shell=False)
		os.chdir('..')
		print("HDECAY in minimal mode terminated.\n")

		# Read the fermion masses from the fermion mass file
		filenameMasses = "HDECAY" + os.sep + "fermionmasses.dat"
		massFileLines = list(line.rstrip('\n') for line in open(filenameMasses))
		MCOSCalc = float((massFileLines[0].split('='))[1].strip())
		MBOSCalc = float((massFileLines[1].split('='))[1].strip())
		
		# Copy the file name to the HDECAY folder and truncate it at the end
		print("Copying input files into HDECAY folder...")
		filenameIn = "Input" + os.sep + inputFileTemp
		filenameOut = "HDECAY" + os.sep + "hdecay.in"
		if os.path.isfile(filenameOut):
			os.remove(filenameOut)
		fileHandler = open(filenameIn, "r")
		convertedFileHandler = []
		lineCount = 1
		for line in fileHandler:
			# If the parameter type is 2, then the reference renormalization scheme has to be set to zero
			if (lineCount == lineWhereRefScheme):
				refScheme = int((line.split())[2])
				if (refScheme == 0):
					print("\nERROR: REFSCHEM=0 is given as input, but REFSCHEM must be greater than zero.")
					print("2HDECAY will be terminated now (ERROR: REFSCHEM=0 was set).")
					sys.exit()
				if (refScheme != 0) and (paramType == 2):
					convertedFileHandler.append("REFSCHEM = 0\n")
					print("\nWARNING: REFSCHEM (line 59) is given as a non-zero value, but TYPE=2 (line 57, lambdas as input) is set.")
					print("REFSCHEM is overwritten to zero, automatic parameter conversion is deactivated!")
					lineCount += 1
					continue
			# Write the current line in an array
			convertedFileHandler.append(line)
			# Check for the parameter type
			if lineCount == lineWhereParamType:
				paramType = int((line.split())[2])
			# Pick out the values of MW, MZ and alphaAtMZ
			if lineCount == lineWhereMZ:
				massMZ = float((line.split())[2])
			if lineCount == lineWhereMW:
				massMW = float((line.split())[2])
			if lineCount == lineWhereAlphaAtMZ:
				alphaAtMZ = float((line.split())[2])
			lineCount += 1
		fileHandler.close()

		# Check whether the new file format (2HDECAY 1.1.0 and upwards) is used
		isNewFormat = False
		for line in convertedFileHandler:
			if "REFSCHEM" in line:
				isNewFormat = True
		if not isNewFormat:
			print("\nERROR: incompatible file format! The input file format was changed with 2HDECAY 1.1.0 to enable input parameter conversion. Please check your input files for compatibility with the new input file format, cf. the Changelog.md file and the manual of the latest 2HDECAY version.")
			print("2HDECAY will be terminated now (ERROR: incompatible file format for 2HDECAY >= 1.1.0).")
			sys.exit()

		# Write a copy of the file to the output folder, but add GFCALC, MCOSCALC and MBOSCALC with the calculated values
		GFcalc = pi/sqrt(2)*alphaAtMZ/(massMW**2*(1-massMW**2/massMZ**2))
		GFline = "GFCALC   = " + str(GFcalc) + "\n"
		MCOSline = "MCOSCALC = " + str(MCOSCalc) + "\n"
		MBOSline = "MBOSCALC = " + str(MBOSCalc) + "\n"
		lineCount = 1
		convertedFile = ''
		renScaleIsDynamic = '0'
		for line in convertedFileHandler:
			if "OUTSCALE" in line and "MIN" in line:
				renScaleIsDynamic = '1'
			if lineCount == lineWhereGFCalc:
				convertedFile += GFline
			elif lineCount == lineWhereOSMC:
				convertedFile += MCOSline
				convertedFile += MBOSline
				convertedFile += line
			else:
				convertedFile += line
			lineCount += 1
		fileHandler = open(filenameOut, "w")
		fileHandler.write(convertedFile)
		fileHandler.close()
		print("... done.\n")
		
		# Calculate the electroweak corrections
		print("Calculating electroweak corrections...\n")
		prompt = ['./electroweakCorrections', '0', '0', '0', '1', 'HDECAY' + os.sep + 'hdecay.in', 'hdecay.in', renScaleIsDynamic]
		subprocess.call(prompt, stdin=None, stdout=None, stderr=None, shell=False)
		print("Calculation of electroweak corrections done.\n")

		# Replace the newline character in each file with a proper newline
		print("Postprocessing temporary input file...")
		fileHandler = open(filenameOut, "r")
		convertedFile = ''
		lineCount = 1
		for line in fileHandler:
			# Convert the literal newlines to actual ones and remove the leading whitespace from the Fortran output
			if (lineCount == lineToInsert):
				lineToReplace = line.replace('\\n', '\n')[1:]
			else:
				lineToReplace = line.replace('\\n', '\n')
			convertedFile += lineToReplace
			# print(lineToReplace)
			lineCount += 1
		fileHandler.close()

		# Store the results file in the correct directory
		if os.path.isfile(filenameOut):
			os.remove(filenameOut)
		fileHandler = open(filenameOut, "w")
		fileHandler.write(convertedFile)
		fileHandler.close()
		print("... done.\n")

		# Start HDECAY in the normal (non-minimal) configuration
		print("Starting HDECAY in standard mode...")
		os.chdir('HDECAY')
		prompt = ['./run']
		subprocess.call(prompt, stdin=None, stdout=None, stderr=None, shell=False)
		os.chdir('..')
		print("HDECAY in standard mode terminated.\n")

		# Copy the output files to the results folder
		print("Copying input files into output folder...")
		filenameIn = "HDECAY" + os.sep + "slha.out"
		filenameOut = "Results" + os.sep + inputFileTemp.replace('.in', '_BR.out')
		copyfile(filenameIn, filenameOut)
		filenameIn = "HDECAY" + os.sep + "ewpartialwidth.out"
		filenameOut = "Results" + os.sep + inputFileTemp.replace('.in', '_EW.out')
		copyfile(filenameIn, filenameOut)
		print("... done.\n")

		# Cleaning 
		if os.path.isfile("HDECAY" + os.sep + "fermionmasses.dat"):
			os.remove("HDECAY" + os.sep + "fermionmasses.dat")
		if os.path.isfile("HDECAY" + os.sep + "alphaandbeta.dat"):
			os.remove("HDECAY" + os.sep + "alphaandbeta.dat")

		print("Corrections for input file " + inputFileTemp + " done.\n")

	# End of program is reached
	print("All calculations finished. Thanks for using 2HDECAY!")

	# Print the end screen 
	print('''
When you use this program please cite:
	The manual for 2HDECAY:
	 - M. Krause, M. Muhlleitner, M. Spira, Comp. Phys. Commun. 246 (2020) 106852, arXiv:1810.00768 (hep-ph)
	The manuals for HDECAY:
	 - A. Djouadi, J. Kalinowski, M. Spira, Comp. Phys. Commun. 108 (1998) 56, hep-ph/9704448
	 - A. Djouadi, J. Kalinowski, M. Muhlleitner, M. Spira, arXiv:1801.09506 (hep-ph)
	The papers on the EW correction to the 2HDM decays:
	 - M. Krause, R. Lorenz, M. Muhlleitner, R. Santos, H. Ziesche, JHEP 1609 (2016) 143, arXiv:1605.04853 (hep-ph)
	 - M. Krause, M. Muhlleitner, R. Santos, H. Ziesche, Phys.Rev. D95 (2017) no.7, 075019, arXiv:1609.04185 (hep-ph)
	 - A. Denner, S. Dittmaier, J.-N. Lang, J. High Energ. Phys. (2018) 2018: 104, arXiv:1808.03466 (hep-ph)
	The publication of LoopTools:
	 - T. Hahn, M. Perez-Victoria, Comp. Phys. Commun. 118 (1999) 153-165, hep-ph/9807565

2HDECAY is released under GNU General Public License (GNU GPL-3.0-or-later). This program is free software: 
you can redistribute it and/or modify it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
for more details.

You have received a copy (LICENSE.md) of the GNU General Public License along with this program in the 2HDECAY
root directory.

Copyright 2018-2019, Marcel Krause, Milada Margarete Muehlleitner and Michael Spira.
	''')

	sys.exit()
