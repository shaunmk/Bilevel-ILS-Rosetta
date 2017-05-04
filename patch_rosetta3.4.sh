#!/bin/bash

### TODO: this whole system could be replaced with one massive diff and patch operation, e.g.
# diff -Naur --exclude=".DS_Store" --exclude="*.pyc" rosetta_source/src rosetta3.4/rosetta_new_sampler/src/
# see man page for patch.
# currently the difficulty is the presence of a number of binary files in the src directory (mostly .pyc ?).

set -e
#set -x

if [ -z $1 ]; then
    echo "Usage: ./patch_rosetta3.4.sh /path/to/rosetta_source" >&2
    exit 1
fi

# SRC_DIR_TO_PATCH is the path to the rosetta3.4 source directory, which will usually end in "rosetta_source"
SRC_DIR_TO_PATCH=$1 # ${HOME}/rosetta3.4_testpatch/rosetta_source
DIR_WITH_PATCHFILES=${PWD}

# check that the supplied path is correct; all the files listed in list_of_patches.txt should be present and have nonzero size.
# This can be obviated by some rigorous check for which version of Rosetta is in use...
echo "***Checking supplied path..."

for patchf in $(cat ${DIR_WITH_PATCHFILES}/list_of_patches.txt)
do
    fname_to_patch=$(basename $patchf .patch)
    dest_dir=$(dirname $patchf)

    target=${SRC_DIR_TO_PATCH}/${dest_dir}/$fname_to_patch

    if [ ! -s $target ]; then
	echo "Error: file $target was not found. Patch cannot be applied." >&2
        echo "  This script will NOT work with versions of Rosetta other than version 3.4." >&2
        echo "  If you are using Rosetta 3.4, there is most likely an error with the path supplied to the script." >&2
	echo "  This will typically end in \"rosetta_source\"." >&2
        exit 2
    fi
done

echo "***Path looks OK."

# copy files that need to be copied. This operation is harmless, even if the subsequent step is unsuccessful.

echo "***Copying new files..."

for copyf in $(cat ${DIR_WITH_PATCHFILES}/filelist_to_copy_not_patch.txt)
do
    dest_dir=$(dirname $copyf)
    target=${SRC_DIR_TO_PATCH}/$dest_dir
    mkdir -p $target
    cp ${DIR_WITH_PATCHFILES}/$copyf $target
done

# patch files that need to be patched. This is a destructive (but reversible!) operation.
# TODO: implement patch-reversal system in case only some files were patched. patch -R

# the same checking as above is done here, because it was done here first.
echo "***Patching existing files..."

for patchf in $(cat ${DIR_WITH_PATCHFILES}/list_of_patches.txt)
do
    fname_to_patch=$(basename $patchf .patch)
    dest_dir=$(dirname $patchf)

    target=${SRC_DIR_TO_PATCH}/${dest_dir}/$fname_to_patch

    if [ -s $target ]; then
	patch $target ${DIR_WITH_PATCHFILES}/$patchf
    else
	echo "Error: file $target was not found or is empty. Patch cannot be applied." >&2
	echo "  This script will NOT work with versions of Rosetta other than version 3.4." >&2
	echo "  If you are using version 3.4, this error really should never be triggered at this point. Contact author." >&2
	exit 3
    fi
done

echo "***Updating Rosetta options list..."

cd $SRC_DIR_TO_PATCH
./update_options.sh

echo "***Patching succeeded!"
echo "    Before compiling, be sure to add your compiler version to tools/build/options.settings if it is not already there."
