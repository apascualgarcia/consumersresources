#!/bin/bash
#!/bin/bash
FOLDER=run_2
SN=2
PATTERN=_opt_mat_intra_specific_syntrophy=allowed_seed_number=2_gamma0=1_alpha0=1

ls $FOLDER/*$PATTERN

#Iterate the loop to read all text files
for value in `ls $FOLDER/*$PATTERN`;
do
    #Read the basename of the file
    filename=`basename $value $PATTERN`
    #Rename all files to doc files
    mv $value $FOLDER/$filename;
done
