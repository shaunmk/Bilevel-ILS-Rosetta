#!/bin/bash

options_file=options.txt
##generates an ${options_file} file with default options for a given protein PDB ID
##this will overwrite any existing file with the same name

if [ -z $1 ]; then 
    echo "${0}: No args; supply protein/target name."
    exit 0
fi

prot=$1

nstruct=1 # Currently, this needs to be set to 1 for Bilevel/ILS protocols.

#if [ ! -z $2 ]; then
#    i=$2
#    if [[ $i =~ ^[0-9]+$ ]]; then
#	nstruct=$i
#    else
#	echo "WARNING: nstruct value was set to ${i} ! Deafulting to nstruct=1."
#    fi
#else
#    echo "WARNING: nstruct value was not supplied! Deafulting to nstruct=1."
#fi

#if [ -z $3 ]; then
#    seed=$RANDOM
#    echo "WARNING: No random seed supplied. Using ${seed}."
#else
#    seed=$3
#fi

# You need to redefine these
fraglib_dir=${HOME}/new_fraglibs
tgt_suffix=_nohoms
database_dir=${HOME}/rosetta3.4/rosetta_database

# to fix the random number seed you'll need these two options
#echo -run:constant_seed >> ${options_file}
#echo -jran $seed >> ${options_file}

# mandatory options

echo -database $database_dir > $options_file
echo -nstruct ${nstruct} >> ${options_file}
echo -frags:annotate >> ${options_file}
echo -abinitio:skip_convergence_check >> ${options_file}
echo -out:pdb >> ${options_file}

# This is the ``production'' setting of increase_cycles for the bilevel/ILS protocols, and was used for all data generated.
# Setting this to 100 uses the same budget of scorefunction evaluations as 10 runs of Rosetta with increase_cycles 10.
echo -abinitio:increase_cycles 100 >> ${options_file}

# This option applies a single round of FastRelax to each structure generated in stage4.
#  A better approach is to run multiple rounds of the procedure separately, using Rosetta's relaxx application.
# echo -abinitio:relax >> ${options_file}

# New options. See README or test/options.txt for explanations
# There are additional options not listed here.
echo -abinitio:SMK_generate_trajectory false >> ${options_file}
echo -abinitio:SMK_var_Linsert_length false >> ${options_file}
echo -abinitio:SMK_SS_dependent_LF false >> ${options_file}

#score function tweaks
#echo -abinitio:stage1_patch ${HOME}/rosetta3.4/rosetta_database/scoring/weights/SMK_weights/score0.wts_SMK >> ${options_file}
#echo -abinitio:stage2_patch ${HOME}/rosetta3.4/rosetta_database/scoring/weights/SMK_weights/score1.wts_SMK >> ${options_file}
#echo -abinitio:stage3a_patch ${HOME}/rosetta3.4/rosetta_database/scoring/weights/SMK_weights/score2.wts_SMK >> ${options_file}
#echo -abinitio:stage3b_patch ${HOME}/rosetta3.4/rosetta_database/scoring/weights/SMK_weights/score5.wts_SMK >> ${options_file}
#echo -abinitio:stage4_patch ${HOME}/rosetta3.4/rosetta_database/scoring/weights/SMK_weights/score3.wts_SMK >> ${options_file}

# some protein names are 4-letter PDB id only, others are PDB id plus chain id
# you may or may not need this
if [ ${#prot} == 4 ]; then
    n=$(echo "${prot}_")
else
    n=$prot
fi

# compose paths to the input files
echo -in:file:fasta ${fraglib_dir}/${prot}${tgt_suffix}/vf_${n}.fasta  >> ${options_file}
echo -in:file:frag3 ${fraglib_dir}/${prot}${tgt_suffix}/${n}.200.3mers >> ${options_file}
echo -in:file:frag9 ${fraglib_dir}/${prot}${tgt_suffix}/${n}.200.9mers >> ${options_file}
echo -psipred_ss2   ${fraglib_dir}/${prot}${tgt_suffix}/${n}.psipred_ss2 >> ${options_file}

# if a native PDB structure is available (e.g. in benchmarking scenarios) this option enables some additional diagnostics.
# No cheating, I promise! :)
#echo -in:file:native ${fraglib_dir}/${prot}${tgt_suffix}/${n}.pdb >> ${options_file}

# mute all stdout
echo -mute all >> ${options_file}

# mute only specific streams
# echo -out:mute:protocols.moves.MonteCarlo >> ${options_file}
# echo -out:mute:protocols.moves.SimulatedAnneal >> ${options_file}
