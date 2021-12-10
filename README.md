# Bio-FlatScope

The codes for our lensless microscope reconstruction (link to be added), proposed in:

**Nature BME**: "In vivo imaging of vasculature and neural activity with a flat, lensless microscope", Jesse K. Adams<sup>&</sup>, Dong Yan<sup>&</sup>, Jimin Wu<sup>&</sup>, Vivek Boominathan<sup>&</sup>, Sibo Gao, Alex V. Rodriguez, Soonyoung Kim, Jennifer Carns, Rebecca Richards-Kortum, Caleb Kemere, Ashok Veeraraghavan* and Jacob T. Robinson*

<sub><sup>&</sup> Denotes equal contribution. | * Corresponding authors </sub>

## 2D Reconstruction
### Data Files
Two supporting data files are necessary for 2D reconstruction:
  1. A single point spread function (PSF) at a certain depth, saved as a ".mat" file.
  2. Raw capture from Bio-FlatScope, saved as an image (.tiff).
### Included Example
We have provided two examples for Bio-FlatScope 2D reconstruction. One is the captured USAF 1951 Resolution test target and the other is a Convallaria majalis slice. Each raw capture has is corresponding PSF. (Results are shown in Figure 1 in our paper).
### Running the code
The code can be tested using our provided examples. A GPU and Matlab Parrallel Computing Toolbox are needed to speed up the reconstruction. 


## 3D Reconstruction
### Data Files
Two supporting data files are necessary for 3D reconstruction:
  1. A stack of point spread functions (PSFs) at different axial distance, saved as a ".mat" file. Available [here](https://drive.google.com/file/d/1UYPXWlYjghcT7DvZNz0ZURw5mnc63Mzf/view?usp=sharing) (Google Drive)
  2. Raw capture from Bio-FlatScope, saved as a ".mat" file.
### Included Example
We have provided one example for Bio-FlatScope 3D reconstruction. The raw data is from a 3D test sample prepared by suspending fluuorescent beads in PDMS phantom. (Results are shown in Figure 2 in our paper)
### Running the code
The code can be tested using our provided examples. A GPU and Matlab Parrallel Computing Toolbox are needed to speed up the reconstruction. 


## Contact Us
In case of any queries regarding the code, please reach out to [Jimin](mailto:jimin.wu@rice.edu) or [Vivek](mailto:vivekb@rice.edu]).
Other raw and analysed data are available for research purpose from corresponding author upon reasonable request.
