# Stellar-cluster-age-from-color-magnitude-diagram
Code and details for CMD estimation of a cluster's age. Observational data taken at the Wallace Astronomical Observatory. Data provided for obtaining identical results for code-testing by users.

Article: https://www.researchgate.net/publication/357924374_Determining_M52's_cluster_age_from_its_Color-Magnitude_Diagram

Run AnalysisProject.py in order to generate plots sequencially. Save them manually and close the windows (the last plot generated will be the CMD plot with Padova isochrones on it). The script assumes that the Padova isochrones text files, the Python script itself and the FIT files (Darks, Biases, Flat fields and cluster images in i' and r' filters) are all in the same directory. The dataset used for generating the article's data is available at: https://www.dropbox.com/s/wth8pes3i322yfv/DATA.zip?dl=0

![alt text](https://github.com/codr-oneci/Stellar-cluster-age-from-color-magnitude-diagram/blob/main/CMD_Isochrones_M52.png)
![alt text](https://github.com/codr-oneci/Stellar-cluster-age-from-color-magnitude-diagram/blob/main/CMD_Isochrones_M52_detail.png)


Library dependencies: numpy, sep, astropy, matplotlib, isochrones. The cluster of choice, M52 is young, having a metalicity of [Fe/H]=-0.11 (slightly smaller than that of the Sun, by definition 0.0). Relevant Padova (theoretical in the ugriz photometric system) isochrones at http://pleiadi.pd.astro.it/.
