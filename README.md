# Bilevel-ILS-Rosetta
_Note: this page is a work in progress._

Patch for the AbinitioRelax application in the Rosetta macromolecular modelling suite. Implements bilevel optimisation and Iterated Local Search (ILS) for improved conformational sampling.

# Requirements
 - Rosetta v3.4 source code (other versions will not work). This requires a license which is freely available to non-commercial users.
 - C++ compiler (tested with GCC 4.8.2 on Linux and llvm-gcc-4.2.1 on OSX)
 - Python 2 (for running SCons; tested with v2.7.5 on Linux and v2.7.12 on OSX)
 - patch (tested with GNU patch v2.7.1 on Linux and patch v2.5.8 on OSX)
 - Linux (preferred) or Mac OSX (see Note for Mac users below)

# IMPORTANT
This patch modifies the source code for the AbinitioRelax application, as well as some core Rosetta classes. I have not tested and cannot guarantee functionality of other applications in the Rosetta suite after the patch has been applied. I recommend that you keep the patched source tree separate from any other installations of Rosetta that you might have.

# Install
1. Obtain Rosetta source distribution and database (version 3.4 ONLY)

	You need to obtain a Rosetta licence from https://www.rosettacommons.org/software and follow the instructions to download Rosetta. Make sure to download version 3.4, as other versions will not work. You will need at least the two downloads marked "required" (Source code and database).

2. Make sure you can compile and run vanilla Rosetta v3.4.

	This step is to exclude any possible errors that might arise due to problems with vanilla Rosetta in combination with your specific system setup. Rosetta supplies a version of SCons, a build system used for compilation.
	```sh
	$ cd /path/to/rosetta_source
	$ ./external/scons-local/scons.py bin mode=release
	```
	On multi-core systems you can supply e.g. ``` -j24 ``` to SCons to use multiple threads (24 in this example).
	In case of any problems with this step, please see the build documentation at https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation.
	Note however, that some paths etc. might be different to those on your systems, as that page refers to more recent versions of Rosetta.

3. Apply the patch and recompile.
	1. Edit patch_rosetta3.4.sh.
	
		Change the value of the ```SRC_DIR_TO_PATCH``` variable to contain the full path to the ```rosetta_source``` directory you downloaded in step 1.
	
	2. Run patch
		```sh
		$ cd /path/to/Bilevel-ILS-Rosetta
		$ ./patch_rosetta3.4.sh
		```
	3. Compile patched version of Rosetta
		```sh
		$ cd /path/to/rosetta_source
		$ ./external/scons-local/scons.py bin mode=release
		```
4. Test

	Sample input and output data are available in ```test/```. These were generated on a Linux system.
	
	Modify ```options.txt``` to correctly indicate the paths to your patched ```rosetta_source``` and ```rosetta_database``` directories.

### Note for Mac users
I have only tested this on OSX versions <= 10.9.5. I have found that I need to use rather old compilers (llvm-gcc-4.2). I did try with Apple's newer Clang compilers (v6.0), but this did not work. I suspect it might be down to adjusting compiler options - you may be able to get ideas from the configuration files for newer versions of Rosetta, which do work fine on up-to-date MacOS systems. This is just something I never spent much time on, but please let me know if you are able to get it to work!
