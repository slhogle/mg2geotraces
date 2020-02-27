<a href="http://www.geotraces.org/"><img src="http://www.geotraces.org/images/banners/GEOTRACES_test.jpg" title="GEOTRACES" alt="GEOTRACES"></a>

# INTEGRATE SIMONS METAGENOMES WITH GEOTRACES IDP17 

Workflow for matching geotraces IDP17 samples to marine metagenomes published [here](https://www.nature.com/articles/sdata2018176)

Also calculates the Deep Chloropyll Maximum depth at each geographic location and determines whether one of the metagenome samples overlaps with the DCM.

## ATTRIBUTION

If any of this code is helpful for you please consider citing this repository. [![DOI](https://zenodo.org/badge/243473977.svg)](https://zenodo.org/badge/latestdoi/243473977)


If you use any of the data mentioned or described in this repository it is ESSENTIAL!! that you properly cite the data sources. These include (but are not limited to):

Biller, S. J., P. M. Berube, K. Dooley, and others. 2018. Marine microbial metagenomes sampled across space and time. [Scientific Data 5: 180176.](https://www.nature.com/articles/sdata2018176)

Schlitzer, R., R. F. Anderson, E. M. Dodas, and others. 2018. The GEOTRACES Intermediate Data Product 2017. [Chem. Geol. 493: 210–223.](https://doi.org/10.1016/j.chemgeo.2018.05.040)

The license for GEOTRACES data does not allow it to be distributed outside of the original. You are encouraged to also cite specific papers describing individual measurement types and/or reach out to those folks to include them as collaborators.

> GEOTRACES Intermediate Data Product (IDP) Download Agreement
> The GEOTRACES programme is keen to ensure that the very significant effort and expertise involved in making trace-element and isotope measurements is acknowledged as fully as possibly in subsequent publications.
>
>Users of the GEOTRACES Intermediate Data Product are expected to abide to the following rules regarding citation
>
>To the greatest extent possible, please cite all relevant publications from researchers that made the measurements you use in your work. Details of publications that should be cited are provided point-by-point in the IDP dataset (in the ODV and ASCII versions) and will be updated on the online database as new papers are published. Where your research deals particularly with data measured by a single group of data originators, you are invited to please contact that group to discuss your work prior to publication and jointly consider the synergy and mutual benefit of co-authorship where appropriate.
>
>Where other constraints prevent citation of all relevant publications, for instance where there is a journal limitation on the maximum number of publications that can be cited, or if the dataset is only used in a minor supportive way, please cite the data compilation itself (as below). In such cases, also please cite any original individual papers that you rely on particularly heavily for your research and interpretations.
>
>Where using data from the IDP2017 product in your publication, please cite the data compilation as: Schlitzer, R., Anderson, R. F., Masferrer Dodas, E, et al., The GEOTRACES Intermediate Data Product 2017, Chem. Geol. (2018), https://doi.org/10.1016/j.chemgeo.2018.05.040 .
>
>Where using data from the IDP2014 product in your publication, please cite the data compilation as: Mawji, E., et al., The GEOTRACES Intermediate Data Product 2014, Mar. Chem. (2015), http://dx.doi.org/10.1016/j.marchem.2015.04.005 .
>
>Users of the GEOTRACES Intermediate Data Product shall not distribute downloaded data to any third party.

## WORKFLOW
To match internal sample identifiers (in form `S0031`) to actual NCBI sequence read archive identifiers you need the files found [here](https://static-content.springer.com/esm/art%3A10.1038%2Fsdata.2018.176/MediaObjects/41597_2018_BFsdata2018176_MOESM325_ESM.zip)

### SECTION GA02
[Match metagenomes to bottle data](bottle2mg/bin/geotraces_GA02_bottle_tidier.md)

[Match metagenomes to the DCM](CTD2mg/bin/geotraces_GA02_CTD_tidier.md)

### SECTION GA03
[Match metagenomes to bottle data](bottle2mg/bin/geotraces_GA03_bottle_tidier.md)

[Match metagenomes to the DCM](CTD2mg/bin/geotraces_GA03_CTD_tidier.md)

### SECTION GA10
[Match metagenomes to bottle data](bottle2mg/bin/geotraces_GA10_bottle_tidier.md)

[Match metagenomes to the DCM](CTD2mg/bin/geotraces_GA10_CTD_tidier.md)

## SECTION GP13
[Match metagenomes to bottle data](bottle2mg/bin/geotraces_GP13_bottle_tidier.md)

[Match metagenomes to the DCM](CTD2mg/bin/geotraces_GP13_CTD_tidier.md)