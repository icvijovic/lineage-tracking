README file

The use of this text file is to provide the user with basic usage of the scripts
required to reproduce the analysis in SI section 2 of Nguyen Ba, Cvijovic and 
Rojas Echenique et al., 2019.

A general description of the pipeline and of its usage can be found in pipeline.txt.

A complete pipeline reproducing the analysis for the YPD population can be found in
example_pipeline.txt. Configuration values are the directory name with the reads 
corresponding to the YPD epochs (here named as YPD_epoch_N).

The pipeline, or any intermediate scripts, can be adapted for other purposes.

OUTPUT

The output of the pipeline reproduces the barcode frequencies file for the YPD 
population found in /data/barcode_frequencies

REQUIREMENTS

Python v3 with the regex library (https://pypi.org/project/regex/).
GCC v6 or above.
Perl v5

CONFIGURATION VALUES

Several of the scripts require specific landing pad configurations or multiplexing
indexes used for the library amplification. The scripts can be modified to accept 
other values.