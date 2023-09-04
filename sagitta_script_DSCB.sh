mkdir fits_std_shorter_range_medium_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,      data_histos_std_shorter_range_medium.root

cd fits_std_shorter_range_medium_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_std_shorter_range_medium_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_std_shorter_range_medium_DSCB
rm -r fits_std_shorter_range_medium_DSCB

###################################

mkdir fits_mc_std_shorter_range_medium_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,      mc_histos_std_shorter_range_medium.root

cd fits_mc_std_shorter_range_medium_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_std_shorter_range_medium_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_mc_std_shorter_range_medium_DSCB
rm -r fits_mc_std_shorter_range_medium_DSCB

###################################

mkdir fits_cvh_shorter_range_medium_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,      data_histos_cvh_shorter_range_medium.root

cd fits_cvh_shorter_range_medium_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_cvh_shorter_range_medium_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_cvh_shorter_range_medium_DSCB
rm -r fits_cvh_shorter_range_medium_DSCB

###################################

mkdir fits_mc_cvh_shorter_range_medium_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,      mc_histos_cvh_shorter_range_medium.root

cd fits_mc_cvh_shorter_range_medium_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_cvh_shorter_range_medium_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_mc_cvh_shorter_range_medium_DSCB
rm -r fits_mc_cvh_shorter_range_medium_DSCB

#####################################

mkdir fits_std_shorter_range_tight_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,      data_histos_std_shorter_range_tight.root

cd fits_std_shorter_range_tight_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_std_shorter_range_tight_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_std_shorter_range_tight_DSCB
rm -r fits_std_shorter_range_tight_DSCB

#######################################

mkdir fits_mc_std_shorter_range_tight_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,      mc_histos_std_shorter_range_tight.root

cd fits_mc_std_shorter_range_tight_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_std_shorter_range_tight_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_mc_std_shorter_range_tight_DSCB
rm -r fits_mc_std_shorter_range_tight_DSCB

#########################################

mkdir fits_cvh_shorter_range_tight_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,     data_histos_cvh_shorter_range_tight.root

cd fits_cvh_shorter_range_tight_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_cvh_shorter_range_tight_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_cvh_shorter_range_tight_DSCB
rm -r fits_cvh_shorter_range_tight_DSCB

#########################################

mkdir fits_mc_cvh_shorter_range_tight_DSCB

#                        sigma_start, low, high, alphaL_start, low, high, nL_start, low, high, alphaR_start, low, high, nR_start, low, high, attempt_no
python sagitta_fitter.py       3.0,    0.2,  5,      1.2,      0.9,  5,     3,    1.7,  20,       1.2,      0.9,  5,      3,    1.7,  20,       11,    mc_histos_cvh_shorter_range_tight.root

cd fits_mc_cvh_shorter_range_tight_DSCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_cvh_shorter_range_tight_DSCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/DSCB/fit_histos_mc_cvh_shorter_range_tight_DSCB
rm -r fits_mc_cvh_shorter_range_tight_DSCB

