mkdir fits_std_shorter_range_medium_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, data_histos_std_shorter_range_medium.root

cd fits_std_shorter_range_medium_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_std_shorter_range_medium_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_std_shorter_range_medium_Voigtian
rm -r fits_std_shorter_range_medium_Voigtian

###################################

mkdir fits_mc_std_shorter_range_medium_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, mc_histos_std_shorter_range_medium.root 

cd fits_mc_std_shorter_range_medium_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_std_shorter_range_medium_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_mc_std_shorter_range_medium_Voigtian
rm -r fits_mc_std_shorter_range_medium_Voigtian

###################################

mkdir fits_cvh_shorter_range_medium_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, data_histos_cvh_shorter_range_medium.root

cd fits_cvh_shorter_range_medium_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_cvh_shorter_range_medium_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_cvh_shorter_range_medium_Voigtian
rm -r fits_cvh_shorter_range_medium_Voigtian

###################################

mkdir fits_mc_cvh_shorter_range_medium_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, mc_histos_cvh_shorter_range_medium.root

cd fits_mc_cvh_shorter_range_medium_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_cvh_shorter_range_medium_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_mc_cvh_shorter_range_medium_Voigtian
rm -r fits_mc_cvh_shorter_range_medium_Voigtian

#####################################

mkdir fits_std_shorter_range_tight_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, data_histos_std_shorter_range_tight.root

cd fits_std_shorter_range_tight_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_std_shorter_range_tight_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_std_shorter_range_tight_Voigtian
rm -r fits_std_shorter_range_tight_Voigtian

#######################################

mkdir fits_mc_std_shorter_range_tight_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, mc_histos_std_shorter_range_tight.root

cd fits_mc_std_shorter_range_tight_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_std_shorter_range_tight_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_mc_std_shorter_range_tight_Voigtian
rm -r fits_mc_std_shorter_range_tight_Voigtian

#########################################

mkdir fits_cvh_shorter_range_tight_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, data_histos_cvh_shorter_range_tight.root

cd fits_cvh_shorter_range_tight_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_cvh_shorter_range_tight_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_cvh_shorter_range_tight_Voigtian
rm -r fits_cvh_shorter_range_tight_Voigtian

#########################################

mkdir fits_mc_cvh_shorter_range_tight_Voigtian

python sagitta_fitter.py     2, 0.0, 10,       2, 0.0, 10,      -0.35 -0.45, -0.25,       -0.20, -0.30, -0.10,        -0.009, -0.4, 0,        4, mc_histos_cvh_shorter_range_tight.root

cd fits_mc_cvh_shorter_range_tight_Voigtian

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_cvh_shorter_range_tight_Voigtian/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/Voigtian/fit_histos_mc_cvh_shorter_range_tight_Voigtian
rm -r fits_mc_cvh_shorter_range_tight_Voigtian
