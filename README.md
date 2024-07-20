# biofilms_on_corrugated_surfaces
This repository contains Jupyter and Mathematica notebooks, Python scripts, and data necessary to recreate all plots from the paper 

"Substrate geometry affects population dynamics in a bacterial biofilm", Witold Postek, Klaudia Staskiewicz, Elin Lilja, and Bartlomiej Waclaw, 
PNAS 121 (17), e2315361121 (2024).
https://www.pnas.org/doi/abs/10.1073/pnas.2315361121

To re-create all figures from the main text and Supplementary Information, (i) download/clone the repository to a local folder, (ii) run all Mathematica and Jupyter notebooks. 
To create only a single chosen figure, run the corresponding notebook only. It is also possible to generate only a single panel from multi-panel figures by running only the relevant part of the notebook.

Some (but not all) notebooks export figures into the data folders; the figures are saved as either PNG or SVG.

Here is a list of all notebooks and the corresponding figures they generate:

Notebook name  | Notebook type | Figures it creates
------------- | ------------- | -------------
Figs_1_S2_S3ABG_S4B.ipynb | Jupyter notebook | Figure 1, SI Figures S2, S3 (panels A,B,G), S4 (panel B)
Figs_2_S13.ipynb | Jupyter notebook | Figure 2, SI Figure S13
Fig_3_S3EF.nb | Wolfram Mathematica | Figure 3, SI Figure S3 panels E,F
Fig_4ABF.nb  | Wolfram Mathematica | Figure 4 panels A,B,F
Figs_4CE_S8_S9.ipynb | Jupyter notebook | Figure 4 panels C,E, SI Figures S8, S9
SI_Fig_S3CD.ipynb | Jupyter notebook | SI Figure S3 panels C,D
SI_Fig_S5.ipynb | Jupyter notebook | SI Figure S5
SI_Fig_S6.ipynb | Jupyter notebook | SI Figure S6
SI_Fig_S7.nb | Wolfram Mathematica | SI Figure S7
SI_Fig_S10AB.ipynb | Jupyter notebook | SI Figures S10 panels A,B
SI_Fig_S11.nb | Wolfram Mathematica | SI Figure S11

The folders included in this repository contain the following:

Folder name | Data type
------------- | -------------
data_Fig_1 | sector boundaries (.csv), heterozygocity (.npy) 
data_Fig_2 | velocity profiles, images of wells (.npb)
data_Fig_3 | computer simulations, no selective advantage, different shapes of wells (.dat)
data_Fig_4 | computer simulations of the experiment from Fig. 4 (.dat), resistant fraction (.txt), sector size data (.csv)
data_S3CD | positions of bacteria at inoculation, flat-bottom wells (.npy)
data_S3EF | computer simulations, different initial densities of bacteria (.dat)
data_S7A | plate reader data (.dat)
data_S7B | plate reader data (.dat)
python_scripts" | Python scripts that have been used to process raw image data (not included here) into data sets included in the folders.

Supplementary Videos that accompany the paper can be found here:

[SI_Video_1.mp4](https://drive.google.com/file/d/14iyic3otvFJQu2msvVp4BJx5odrc9IEq/view?usp=drive_link) - biofilm growth in selected wells (as in Fig. 1) from inoculation (t=0) until t=144 h.

[SI_Video_2.mp4](https://drive.google.com/file/d/16Ci4M7xqwNQUY3_ruPsHmT-WE9_3XJOf/view?usp=drive_link) - biofilm growth in 4 wells of different type: flat, (T,A)=(10,1.7), (T,A)=(20,5.1), and (T,A)=(50,8.7), from inoculation (t=0) until t=182 h.

[SI_Video_3.mp4](https://drive.google.com/file/d/1HvDA5ihCHf6GBohLbYQ-m3IseaWurSmb/view?usp=drive_link) - simulated biofilm in an 80x70 um well, for the same parameters as simulations presented in Fig. 4.  (raw data in test_80x80_normal)

[SI_Video_4.mp4](https://drive.google.com/file/d/1Ubm9To1DWivnUAsXomiZgopGF9U2I8Wb/view?usp=drive_link) - simulated biofilm in an 80x30 um well, for the same parameters as SI_Video_3.mp4. (raw data in test_80x30)

[SI_Video_5.mp4](https://drive.google.com/file/d/1G2bwqpoWgxEBNsyiL8_jH27fp4jQ6j4I/view?usp=drive_link) - simulated biofilm in an 80x30 um well, for the same parameters as SI_Video_3.mp4 with the exception of the friction coefficient being 20x larger. (raw data in test_80x30_zeta_20x)

