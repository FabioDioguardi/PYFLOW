PYFLOW v 2.2
 
Copyright (C) 2020 Fabio Dioguardi[1]
[1] British Geological Survey, The Lyell Centre, Edinburgh, United Kingdom. fabiod@bgs.ac.uk


The package comes with all the routines presented in the previous section, each one in a separate .f90 file. The user needs to compile the FORTRAN90 files and build the executable. In order to simplify this operation, a script (named Makefile) is also provided. The script can be invoked with the freeware Gnu Make software. The user should only run the Make program in the folder in which all the source files and the script are stored by typing make. The command make clean deletes some files created during the compilation: .mod and .o.

- Linux

In Linux operating systems Make should be installed by default, otherwise the user can download and install the program with the package manager specific of the OS or by typing the proper command on the command shell (e.g. apt-get install for Ubuntu, yum for Fedora, etc.). The command which make gives information on whether and where Make is installed.
The Makefile is written assuming that Gfortran compiler is used. If this is not the case, the user can edit the second line of the Makefile by replacing "gfortran" with the proper command invoking the desired compiler. 

- Windows and Mac OS

For Windows and Mac operating systems the user can find the make executable on Internet (e.g. Make for Windows). For these OSs only the executable is available, which has to be placed in the same folder where the Makefile script and the source code are. The other possibility is to work with a Linux emulator (e.g. Cygwin for Windows).

