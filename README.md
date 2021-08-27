# README #


### *Summary:* ###

Jupyter Notebooks & data for recreating the main analyses in [Dimitrova et al., (2021)](https://www.sciencedirect.com/science/article/pii/S1053811921007618) Preterm birth alters the development of cortical microstructure and morphology at term-equivalent age. For further details on pre-processing and analysis, please refer to the paper. 


MRI data includes surface measures for 259 term-born and 76 preterm infants scanned at term-equivalent age for the [developing Human Connectome Project](http://www.developingconnectome.org/), preprocessed using the [dHCP structural](https://github.com/BioMedIA/dhcp-structural-pipeline) [(Makropoulos et al., 2018)](https://pubmed.ncbi.nlm.nih.gov/29409960/) & dHCP diffusion [(Christiaens et al., 2021)](https://www.sciencedirect.com/science/article/pii/S1053811920309228) pipelines. 


*The 3rd release of the dHCP imaging data is now publicly available [here](http://www.developingconnectome.org/data-release/third-data-release/).*

-------

### *Dependencies:* ###

* python => 3.7 
* [Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install.html) 
* Used packages include: `matplotlib` , `nibabel` , `numpy` , `pandas` , `scikit-learn` , `scipy` , `seaborn` , `skopt` , `statsmodels` , [GPy](https://github.com/SheffieldML/GPy) 
* [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench), [Toblerone](https://toblerone.readthedocs.io/en/latest/)


-------

* `./data` - csv files containing parcel-averaged surface measures (FA, MD, ODI, fICVf, cortical thickness, curvature, sulcation & surface area) for all term & preterm infants. 
* `./random_parcels` - left/right Voronoi parcellation surface files (n=143 per hemisphere).
* `./pvc` - partial volume estimation & correction


-------

### *References:* ###

Christiaens, D. et al., ["Scattered slice SHARD reconstruction for motion correction in multi-shell diffusion MRI"](https://www.sciencedirect.com/science/article/pii/S1053811920309228), Neuroimage 225 (2021): 117437.

Dimitrova, R. et al. ["Preterm birth alters the development of cortical microstructure and morphology at term-equivalent age"](https://www.sciencedirect.com/science/article/pii/S1053811921007618), NeuroImage 243 (2021): 118488.

Makropoulos, A. et al. ["The developing human connectome project: A minimal processing pipeline for neonatal cortical surface reconstruction"](https://www.sciencedirect.com/science/article/abs/pii/S1053811918300545?via%3Dihub), Neuroimage 173 (2018): 88-112.

------

### *License:* ###

The data distributed with this work are made available under a [CC-BY 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
