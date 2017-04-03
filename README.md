# Bilevel-ILS-Rosetta
_Note: this page is a work in progress._

Patch for the AbinitioRelax application in the Rosetta macromolecular modelling suite. Implements bilevel optimisation and Iterated Local Search (ILS) for improved conformational sampling.

# Requirements
 - Rosetta v3.4 source code (other versions will not work). This requires a license which is freely available to non-commercial users.
 - C++ compiler (tested with GCC 4.8.2 on Linux and llvm-gcc-4.2.1 on OSX)
 - Python 2 (for running SCons; tested with v2.7.5 on Linux and v2.7.12 on OSX)
 - patch (tested with GNU patch v2.7.1 on Linux and patch v2.5.8 on OSX)
 - Linux (preferred) or Mac OSX (see Note for Mac users below)

## IMPORTANT
This patch modifies the source code for the AbinitioRelax application, as well as some core Rosetta classes. I have not tested and cannot guarantee functionality of other applications in the Rosetta suite after the patch has been applied. I recommend that you keep the patched source tree separate from any other installations of Rosetta that you might have.

# Installation instructions
1. Obtain Rosetta source distribution and database (version 3.4)

	You need to obtain a Rosetta licence from https://www.rosettacommons.org/software and follow the instructions to download Rosetta. Make sure to download version 3.4, as other versions will not work. You will need at least the two downloads marked "required" (Source code and database).

2. Apply the patch.   
	The script ```patch_rosetta3.4``` handles the patch process. Supply it with the full path to the ```rosetta_source``` directory you downloaded in step 1.

	```sh
	$ cd /path/to/Bilevel-ILS-Rosetta
	$ ./patch_rosetta3.4.sh /path/to/rosetta_source
	```
	The script will do some simple checks to ensure that the patch can be applied. 
		
3. Compile patched version of Rosetta.   
	Once the patch has succeeded, you will need to compile the patched Rosetta source tree. Rosetta comes with a version of SCons, a Python-based build system.

	```sh
	$ cd /path/to/rosetta_source
	$ ./external/scons-local/scons.py bin mode=release
	```

	Optionally, you can specify the build target as ```bin/AbinitioRelax``` to only compile everything necessary for that application.
	The compilation process can be quite slow, but on multi-core systems you can supply e.g. ``` -j24 ``` to get SCons to use multiple threads (24 in this example). Successful compilation is usually indicated by
	```
	scons: done building targets.
	```

4. Test (TODO)

	Sample input and output data are available in ```test/```. These were generated on a Linux system.
	
	Modify ```options.txt``` to correctly indicate the paths to your patched ```rosetta_source``` and ```rosetta_database``` directories.

### Common issues
By default, SCons will try to use whatever ```gcc``` refers to in your shell. On standard configurations of MacOS for example, this is usually Apple's LLVM compiler. To change the default, edit the file ```/path/to/rosetta_source/tools/build/basic.options```.

Another common error when running SCons is e.g.
```
### KeyError: "Unknown version number 4.8 for compiler 'gcc'"
```
In this case, you will need to let SCons know about this version of GCC. Edit the file ```/path/to/rosetta_source/tools/build/options.settings``` and add the correct version (in this example, "4.8") to the list of GCC versions provided.

### Note for Mac users
I have only tested this on OSX versions <= 10.9.5. I have found that I need to use rather old compilers (llvm-gcc-4.2, available from MacPorts as apple-gcc42). I did try with Apple's newer LLVM compilers (v6.0), but this did not work. I suspect it might be down to adjusting compiler options in ```/path/to/rosetta_source/tools/build/basic.settings``` - you may be able to get ideas from the version of this file in newer releases of Rosetta, which do work fine on up-to-date MacOS systems. All my production work was done on Linux with GCC, so this is just something I never spent much time on. Please let me know if you are able to get it to work!
