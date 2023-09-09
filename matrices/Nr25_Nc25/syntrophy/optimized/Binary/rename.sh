#!/bin/bash
#!/bin/bash
FOLDER=run_3
SN=3

ls $FOLDER/*_opt_mat_intra_specific_syntrophy=allowed_seed_number=$SN

#Iterate the loop to read all text files
for value in `ls $FOLDER/*_opt_mat_intra_specific_syntrophy=allowed_seed_number=$SN`;
do
    #Read the basename of the file
    filename=`basename $value _opt_mat_intra_specific_syntrophy=allowed_seed_number=$SN`
    #Rename all files to doc files
    mv $value $FOLDER/$filename;
done
