#mkdir fits_mc_ideal_shorter_range_medium_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,        mc_histos_ideal_shorter_range_medium.root

#cd fits_mc_ideal_shorter_range_medium_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_mc_ideal_shorter_range_medium_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_mc_ideal_shorter_range_medium_BWCB
#rm -r fits_mc_ideal_shorter_range_medium_BWCB

###################################

#mkdir fits_std_shorter_range_medium_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,        data_histos_std_shorter_range_medium.root

#cd fits_std_shorter_range_medium_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_std_shorter_range_medium_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_std_shorter_range_medium_BWCB
#rm -r fits_std_shorter_range_medium_BWCB

###################################

mkdir fits_mc_std_shorter_range_medium_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,        mc_histos_std_shorter_range_medium.root

cd fits_mc_std_shorter_range_medium_BWCB

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
cd ..

cp /home/users/alexe/workingarea/Sagitta/fits_mc_std_shorter_range_medium_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_mc_std_shorter_range_medium_BWCB
rm -r fits_mc_std_shorter_range_medium_BWCB

###################################

#mkdir fits_cvh_shorter_range_medium_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,       data_histos_cvh_shorter_range_medium.root

#cd fits_cvh_shorter_range_medium_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_cvh_shorter_range_medium_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_cvh_shorter_range_medium_BWCB
#rm -r fits_cvh_shorter_range_medium_BWCB

###################################

#mkdir fits_mc_cvh_shorter_range_medium_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,       mc_histos_cvh_shorter_range_medium.root

#cd fits_mc_cvh_shorter_range_medium_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_mc_cvh_shorter_range_medium_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_mc_cvh_shorter_range_medium_BWCB
#rm -r fits_mc_cvh_shorter_range_medium_BWCB

#####################################

#mkdir fits_std_shorter_range_tight_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,        data_histos_std_shorter_range_tight.root

#cd fits_std_shorter_range_tight_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_std_shorter_range_tight_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_std_shorter_range_tight_BWCB
#rm -r fits_std_shorter_range_tight_BWCB

#######################################

#mkdir fits_mc_std_shorter_range_tight_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1, mc_histos_std_shorter_range_tight.root

#cd fits_mc_std_shorter_range_tight_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_mc_std_shorter_range_tight_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_mc_std_shorter_range_tight_BWCB
#rm -r fits_mc_std_shorter_range_tight_BWCB

#########################################

#mkdir fits_cvh_shorter_range_tight_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,     data_histos_cvh_shorter_range_tight.root

#cd fits_cvh_shorter_range_tight_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_cvh_shorter_range_tight_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_cvh_shorter_range_tight_BWCB
#rm -r fits_cvh_shorter_range_tight_BWCB

#########################################

#mkdir fits_mc_cvh_shorter_range_tight_BWCB

#                         width_start, low, high, sigma_start, low, high, alpha_start, low, high, n_start, low, high, attempt_no
#python sagitta_fitter.py       2, 0.0, 10,           2, 0.0, 10,                 1,  0.1, 10,       10, 3,  15,   1,    mc_histos_cvh_shorter_range_tight.root

#cd fits_mc_cvh_shorter_range_tight_BWCB

#gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_fits_6_6_6_6.pdf -dBATCH sim_fit_region_*_*_*_*.pdf
#cd ..

#cp /home/users/alexe/workingarea/Sagitta/fits_mc_cvh_shorter_range_tight_BWCB/combined_fits_6_6_6_6.pdf /home/users/alexe/workingarea/Sagitta/BWCB/fit_histos_mc_cvh_shorter_range_tight_BWCB
#rm -r fits_mc_cvh_shorter_range_tight_BWCB
