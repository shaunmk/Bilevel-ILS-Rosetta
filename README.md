# Bilevel-ILS-Rosetta

Patch for the AbinitioRelax application in the Rosetta macromolecular modelling suite. Implements bilevel optimisation and Iterated Local Search (ILS) for improved conformational sampling.

# Prerequisites
 - Rosetta v3.4 source code (other versions will not work). This requires a license which is freely available to non-commercial users.
 - C++ compiler (tested with GCC 4.8.2 on Linux and llvm-gcc-4.2.1 on OSX)
 - Python 2 (for running SCons; tested with v2.7.5 on Linux and v2.7.12 on OSX)
 - patch (tested with GNU patch v2.7.1 on Linux and patch v2.5.8 on OSX)
 - Linux (preferred) or Mac OSX (see Note for Mac users below)

# IMPORTANT
This patch modifies the source code for the AbinitioRelax application, as well as some core Rosetta classes. I have not tested and cannot guarantee functionality of other applications in the Rosetta suite after the patch has been applied. I recommend that you keep the patched source tree separate from any other installations of Rosetta that you might have.

The ```rosetta_database/``` directory is not modified in any way, so you can keep a single copy of that.

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

4. Test

	Sample inputs and outputs are available in ```test/```. These were generated on a Linux system.
	
	Modify ```test/inputs/options.txt``` to correctly indicate the paths to your patched ```rosetta_source``` and ```rosetta_database``` directories. Then run the patched executable:
	
	```sh
	$ cd /path/to/Bilevel-ILS-Rosetta/test/inputs
	$ /path/to/rosetta_source/bin/AbinitioRelax.default.linuxgccrelease @options.txt
	```
	
	Note that the extension applied to the executable will depend on the system and compiler used. For example, on a Mac using clang it would be ```AbinitioRelax.default.macosclangrelease```. As far as I'm aware, there is no difference between the ```.default``` and non-```.default``` versions of the executable.

### Other directories
- The ```example``` directory   
	This directory contains an example options file and a helper script that generates it. The settings in this options file are typical of a "production" run, and these settings were used for all results in our paper. The key difference between these settings and those in ```test/inputs/options.txt``` is the setting of the ```increase_cycles``` parameter.

- The ```tools``` directory   
	This directory contains the script used to make the patch files, and a list of all files that were added or modified relative to standard Rosetta. You can use the script if you want to make your own version of our protocols. You may need to modify paths defined therein.   

### Common issues with building
By default, SCons will try to use whatever ```gcc``` refers to in your shell. On standard configurations of MacOS for example, this is usually Apple's LLVM compiler. To change the default, edit the file ```/path/to/rosetta_source/tools/build/basic.options```.

Another common error when running SCons is e.g.
```
### KeyError: "Unknown version number 4.8 for compiler 'gcc'"
```
In this case, you will need to let SCons know about this version of GCC. Edit the file ```/path/to/rosetta_source/tools/build/options.settings``` and add the correct version (in this example, "4.8") to the list of GCC versions provided.

### Note for Mac users
I have only tested this on OSX versions <= 10.9.5, using rather old compilers (llvm-gcc-4.2). I think this version was current at the time Rosetta 3.4 was released, and I simply kept my system that way. The easiest way to get up and running is to use this (or similar) compiler as your system default ```gcc```. If you do not want this, it is still possible to use a compiler other than the system default. I was able to use ```gcc48``` from macports, but this requires editing the build configuration files in ```/path/to/rosetta_source/tools/build```. Please raise an issue and let me know if you are trying to do this, as it is a bit fiddly.

I did try to build with Apple's newer LLVM compilers (v6.0), but this did not work. I suspect it might be down to adjusting compiler options in ```/path/to/rosetta_source/tools/build/basic.settings``` - you may be able to get ideas from the version of this file in newer releases of Rosetta, which do work fine on up-to-date MacOS systems. All my production work was done on Linux with GCC, so this is just something I never spent much time on. Please let me know if you are able to get it to work!

# Running
Follow the approach in Step 4 of the Installation instructions. An example "production" version of the options file is given in the ```example/``` directory. To generate sets of 1000 decoys, we ran 100 independent runs of our protocols with these options, the key being ```-abinitio:increase_cycles 100```. This setup uses the same budget of function evaluations as 1000 runs of regular Rosetta AbinitioRelax with ```-abinitio:increase_cycles 10```.

The fragment sets used in our papers can be found here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1254031.svg)](https://doi.org/10.5281/zenodo.1254031)

# Configuring custom archivers
In addition to the default energy-based archiving strategy, you can choose to use a couple of different archiving strategies for local minima. Two examples are the Stochastic Ranking-based archiver operating on contact maps (SRCM) and the Elitist Random (ER) archiver. For details on these archivers, consult Chapter 5 of my PhD thesis, which is available [here](https://www.research.manchester.ac.uk/portal/en/theses/on-conformational-sampling-in-fragmentbased-protein-structure-prediction(15b113b0-c389-4e8e-ae7c-d399e8a72b03).html). A paper based on this chapter is in preparation.

Once you have applied the bilevel-ILS patch to Rosetta, navigate to `/path/to/rosetta_source/src/protocols/moves/MonteCarlo.cc` and edit the (mis-named) function `MonteCarlo::create_archives_from_cmdline()`, using the examples provided in the comments there (see also Known Isuues below). Some archiver types have additional options, examples of which are also provided. Each archiver is added to a `std::vector < MGFUtils::MGFArchive >`, and has the prototype 

`MGFUtils::MGFArchive(desired_size, max_size, archive_type, tag)`, where 

`desired_size` and `max_size` are `int`s that specify the desired and maximum number of structures in the archive, `archive_type` is one of a few predefined strings that selects the type of archiver, and `tag` is a string that is used to uniquely identify the specific archiver, and is appended to the filenames of the output structures. `tag` must be unique across all archivers. Some archivers can be configured with additional options, examples of these are given. An example is the 'rho' parameter for archives based on stochastic ranking, which can be set e.g. to 0.5 using

`Archives.back().set_SR_PROB(0.5);`

immediately after adding the archive.

You'll notice that a few additional archivers are also implemented; feel free to experiment! Only the energy-based archiver is enabled by default. Once you've made your edits, you will need to recompile the codebase using the instructions above.

# Citing
If you find this repo useful, please cite our [paper in Scientific Reports](https://www.nature.com/articles/s41598-018-31891-8).

#### TODO
- Table with all newly defined command-line options, what they do and their defaults

# Known issues
- The ```nstruct``` parameter must always be set to 1.   

	Setting this parameter to anything greater than 1 will still _work_, but you will overwrite output files from prior iterations of the ```nstruct``` loop. In order to run the application many times to get a larger decoy set, I recommend running each run in a separate working directory.   
	
	The actual number of structures returned by the application is set when configuring the solution archiving strategy(ies) employed. This is done in the (admittedly mis-named) function ```MonteCarlo::create_archives_from_cmdline()``` in ```src/protocols/moves/MonteCarlo.cc```. My plan was eventually to introduce some sort of configuration file for this part of the application, so that users don't need to recompile everytime they want to try out a new configuration of archivers. Please raise an issue if this is something you'd find useful.

- The commandline parameter ```-frags:annotate``` must (currently) always be provided.   

	This option makes the application keep track of the source of the structural info (PDB template and residue numbers) for each residue in the structure, following each fragment insertion. Correctly implementing this was tricky and needed many edits to many files.


