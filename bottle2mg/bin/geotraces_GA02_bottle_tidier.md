GA02 bottle tidier
================

  - [Preprocessing of the discrete bottle file using GNU unix
    tools](#preprocessing-of-the-discrete-bottle-file-using-gnu-unix-tools)
  - [Processing with R](#processing-with-r)
  - [Read metagenomes metadata](#read-metagenomes-metadata)
  - [Read geotraces IDP17 discrete bottle data by cruise
    section](#read-geotraces-idp17-discrete-bottle-data-by-cruise-section)
  - [Trace metal data from the
    bottles/pumps](#trace-metal-data-from-the-bottlespumps)
  - [Getting data tidied](#getting-data-tidied)

This uses the GEOTRACES Intermediate Data Product 2017 (Version 2)

The license for this data does not allow it to be distributed outside of
the original. See here:

> GEOTRACES Intermediate Data Product (IDP) Download Agreement The
> GEOTRACES programme is keen to ensure that the very significant effort
> and expertise involved in making trace-element and isotope
> measurements is acknowledged as fully as possibly in subsequent
> publications.
> 
> Users of the GEOTRACES Intermediate Data Product are expected to abide
> to the following rules regarding citation
> 
> To the greatest extent possible, please cite all relevant publications
> from researchers that made the measurements you use in your work.
> Details of publications that should be cited are provided
> point-by-point in the IDP dataset (in the ODV and ASCII versions) and
> will be updated on the online database as new papers are published.
> Where your research deals particularly with data measured by a single
> group of data originators, you are invited to please contact that
> group to discuss your work prior to publication and jointly consider
> the synergy and mutual benefit of co-authorship where appropriate.
> 
> Where other constraints prevent citation of all relevant publications,
> for instance where there is a journal limitation on the maximum number
> of publications that can be cited, or if the dataset is only used in a
> minor supportive way, please cite the data compilation itself (as
> below). In such cases, also please cite any original individual papers
> that you rely on particularly heavily for your research and
> interpretations.
> 
> Where using data from the IDP2017 product in your publication, please
> cite the data compilation as: Schlitzer, R., Anderson, R. F.,
> Masferrer Dodas, E, et al., The GEOTRACES Intermediate Data Product
> 2017, Chem. Geol. (2018),
> <https://doi.org/10.1016/j.chemgeo.2018.05.040> .
> 
> Where using data from the IDP2014 product in your publication, please
> cite the data compilation as: Mawji, E., et al., The GEOTRACES
> Intermediate Data Product 2014, Mar. Chem. (2015),
> <http://dx.doi.org/10.1016/j.marchem.2015.04.005> .
> 
> Users of the GEOTRACES Intermediate Data Product shall not distribute
> downloaded data to any third party.

To follow along you muust register with the [the British Oceanographic
Data Centre](https://www.bodc.ac.uk/geotraces/data/idp2017/) and
download the data. We will be using the `Discrete Sample Data` in ODV
format.

Note that there [are some
issues](https://www.bodc.ac.uk/data/documents/nodb/544232/) with some of
the IDP data. Fortunately, none of these issues apply to the samples
corresponding to the metagenomes.

## Preprocessing of the discrete bottle file using GNU unix tools

Split into individual cruise transects. Not sure if this is necessary
but I thought that maybe columns are different between cruise sections.

``` bash
grep "^Cruise" ../../GEOTRACES_IDP2017_v2_Discrete_Sample_Data/GEOTRACES_IDP2017_v2_Discrete_Sample_Data_c2f9b27_3.txt > ../../GEOTRACES_IDP2017_v2_Discrete_Sample_Data/GA02_bottle.tsv
grep "^GA02" ../../GEOTRACES_IDP2017_v2_Discrete_Sample_Data/GEOTRACES_IDP2017_v2_Discrete_Sample_Data_c2f9b27_3.txt >> ../../GEOTRACES_IDP2017_v2_Discrete_Sample_Data/GA02_bottle.tsv
```

## Processing with R

``` r
library(here)
library(tidyverse)
library(janitor)
library(fuzzyjoin)
```

## Read metagenomes metadata

``` r
GA02.mg <- read_tsv(here("tara_simons_timeseries_basic_meta.tsv")) %>% 
  filter(section=="GA02") %>% 
  select(sampleID, bottle_num_BODC, station_num_BODC, lat, lon, depth) %>%
  mutate(station_num_BODC=as.numeric(station_num_BODC)) %>%
  mutate(bottle_num_BODC=as.numeric(bottle_num_BODC)) %>%
  mutate(depthmin=depth-5,
         depthmax=depth+5) %>%
  mutate(depthmin=ifelse(depthmin < 0, 0, depthmin))
```

## Read geotraces IDP17 discrete bottle data by cruise section

``` r
GA02.btl <- read_tsv(here("GEOTRACES_IDP2017_v2_Discrete_Sample_Data", 
              "GA02_bottle.tsv"), col_names = TRUE) %>%
  janitor::clean_names() %>%
  janitor::remove_empty(which = "cols") %>%
  select(-contains("qv_iode"), 
         -contains("infos_"), 
         -contains("standard_dev_"), 
         -qv_odv_sample) %>%
  rename(station_num_BODC=bodc_station, 
         bottle_num_BODC=bodc_bottle_number,
         depth=depth_m) %>%
  select(-type, -yyyy_mm_dd_thh_mm_ss_sss, -bot_depth_m, -operators_cruise_name, -cruise_report, -ship_name, -period, -chief_scientist, -geotraces_scientist, - cast_identifier, -bottle_number, -bottle_flag, -firing_sequence, -infos)
```

## Trace metal data from the bottles/pumps

``` r
GA02.tm <- fuzzy_inner_join(GA02.mg, GA02.btl, 
                      by=c(
                        "station_num_BODC" = "station_num_BODC",
                        "depthmin" = "depth",
                        "depthmax" = "depth"),
                       match_fun = list(`==`, `<=`, `>=`)
                      ) %>%
  janitor::remove_empty(which = "cols") %>%
  mutate(depth_diff=abs(depth.x-depth.y))
```

### some double checks

all mg samples are present

``` r
distinct(GA02.tm, sampleID) %>% nrow()
```

    ## [1] 127

available columns/data

``` r
colnames(GA02.tm)
```

    ##  [1] "sampleID"                        "bottle_num_BODC.x"              
    ##  [3] "station_num_BODC.x"              "lat"                            
    ##  [5] "lon"                             "depth.x"                        
    ##  [7] "depthmin"                        "depthmax"                       
    ##  [9] "cruise"                          "station"                        
    ## [11] "longitude_degrees_east"          "latitude_degrees_north"         
    ## [13] "station_num_BODC.y"              "pressure_dbar"                  
    ## [15] "depth.y"                         "sampling_device"                
    ## [17] "bottle_num_BODC.y"               "ctdtmp_deg_c"                   
    ## [19] "ctdsal"                          "ctdoxy_umol_kg"                 
    ## [21] "phosphate_d_conc_bottle_umol_kg" "silicate_d_conc_bottle_umol_kg" 
    ## [23] "nitrate_d_conc_bottle_umol_kg"   "nitrite_d_conc_bottle_umol_kg"  
    ## [25] "dic_d_conc_bottle_umol_kg"       "doc_d_conc_bottle_umol_kg"      
    ## [27] "al_d_conc_bottle_nmol_kg"        "co_d_conc_bottle_pmol_kg"       
    ## [29] "fe_d_conc_bottle_nmol_kg"        "mn_d_conc_bottle_nmol_kg"       
    ## [31] "mo_d_conc_bottle_nmol_kg"        "ni_d_conc_bottle_nmol_kg"       
    ## [33] "pb_d_conc_bottle_pmol_kg"        "zn_d_conc_bottle_nmol_kg"       
    ## [35] "depth_diff"

## Getting data tidied

This is a pain, but I think it is the only way to really clean this up

``` r
GA02.long.0 <- GA02.tm %>%
  select(sampleID, depth.mg=depth.x, depth.btl=depth.y, depth_diff, pressure_dbar,
         ctdtmp_deg_c:zn_d_conc_bottle_nmol_kg) %>%
  pivot_longer(
    cols = pressure_dbar:zn_d_conc_bottle_nmol_kg,
    names_to = "name", 
    values_to = "value") %>%
  mutate(phase=case_when(         str_detect(name, "ctdtmp_deg_c")  ~ "dissolved",
                                  str_detect(name, "ctdsal")  ~ "dissolved",
                                  str_detect(name, "ctdoxy_umol_kg")  ~ "dissolved",
                                  str_detect(name, "pressure_dbar")  ~ "dissolved",
                                  str_detect(name, "_d_")  ~ "dissolved",
                                  str_detect(name, "_s_")  ~ "soluble",
                                  str_detect(name, "_tp_") ~ "total.particulate")) %>%
  mutate(biogeo_var=case_when(    str_detect(name, "ctdtmp_deg_c")  ~ "temperature",
                                  str_detect(name, "ctdsal")  ~ "salinity",
                                  str_detect(name, "ctdoxy_umol_kg")  ~ "oxygen",
                                  str_detect(name, "pressure_dbar")  ~ "pressure",
                                  str_detect(name, "^phosphate_")  ~ "phosphate",
                                  str_detect(name, "^silicate_")  ~ "silicate",
                                  str_detect(name, "^nitrate_")  ~ "nitrate",
                                  str_detect(name, "^nitrite_")  ~ "nitrite",
                                  str_detect(name, "^doc_") ~ "dissolved.organic.carbon",
                                  str_detect(name, "^dic_") ~ "dissolved.inorganic.carbon",
                                  str_detect(name, "^al_") ~ "aluminum",
                                  str_detect(name, "^co_") ~ "cobalt",
                                  str_detect(name, "^cu_") ~ "copper",
                                  str_detect(name, "^fe_s_|^fe_tp_|^fe_d_") ~ "iron",
                                  str_detect(name, "^fe_ii_d") ~ "ironII",
                                  str_detect(name, "^mn_") ~ "manganese",
                                  str_detect(name, "^mo_") ~ "molybdenum",
                                  str_detect(name, "^ni_") ~ "nickel",
                                  str_detect(name, "^pb_") ~ "lead",
                                  str_detect(name, "^zn_") ~ "zinc",
                                  str_detect(name, "^l1cu_d_conc")   ~ "L1.copper.conc",
                                  str_detect(name, "^l1cu_d_log_k")  ~ "L1.copper.bind",
                                  str_detect(name, "^l1fe_d_conc")   ~ "L1.iron.conc",
                                  str_detect(name, "^l1fe_d_log_k")  ~ "L1.iron.bind",
                                  str_detect(name, "^l2fe_d_conc")   ~ "L2.iron.bind",
                                  str_detect(name, "^l2fe_d_log_k")  ~ "L2.iron.logk")) %>%
  mutate(sample_method=case_when( str_detect(name, "_fish")  ~ "fish",
                                  str_detect(name, "_bottle")  ~ "bottle",
                                  TRUE ~ "bottle")) %>%
  mutate(unit=case_when(          str_detect(name, "_umol_kg")  ~ "umol.kg",
                                  str_detect(name, "_nmol_kg")  ~ "nmol.kg",
                                  str_detect(name, "_pmol_kg")  ~ "pmol.kg",
                                  str_detect(name, "d_log_k")  ~ "logK",
                                  str_detect(name, "_deg_c")  ~ "deg.c",
                                  str_detect(name, "ctdsal")  ~ "PSS.1978",
                                  str_detect(name, "_dbar")  ~ "dbar")) %>%
  mutate(sensitivity=      ifelse(str_detect(name, "_ll_d"), "high", "normal"))
```

Collapse samples with multiple observations and complete missing
observations with `NA`

``` r
GA02.long.1 <- GA02.long.0 %>%
  filter(!is.na(value)) %>% 
  group_by(sampleID, biogeo_var, phase, sensitivity) %>%
  arrange(depth_diff) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(sampleID, biogeo_var, phase) %>%
  arrange(sensitivity) %>%
  slice(1) %>%
  ungroup %>%
  complete(sampleID, nesting(biogeo_var, phase)) %>%
  group_by(biogeo_var) %>%
  fill(unit, .direction="updown") %>% 
  ungroup()
```

Recombine with the metagenome information

``` r
GA02.long.2 <- left_join(GA02.long.1, GA02.mg) %>%
  select(sampleID, station_num_BODC, lat, lon, depth.mg=depth, depth.btl, depth_diff, biogeo_var, unit, phase, sample_method, value)
```

    ## Joining, by = "sampleID"

This is the step you would take to put data into wide format. THis is
what you need for example if you are using the data with statistical
models in R.

``` r
GA02.wide <- GA02.long.1 %>%
  select(-name, -depth.mg, -sensitivity) %>%
  unite(newname, biogeo_var, phase, unit) %>%
  select(-depth.btl, -depth_diff, -sample_method) %>%
  pivot_wider(names_from = "newname", values_from="value") 
```

Finally write out the completed tidied result

``` r
write_tsv(GA02.long.2, here("bottle2mg", "GA02_mg_bottles_long.tsv"))
```

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] fuzzyjoin_0.1.5    janitor_1.2.1      forcats_0.4.0      stringr_1.4.0     
    ##  [5] dplyr_0.8.3        purrr_0.3.3        readr_1.3.1        tidyr_1.0.0       
    ##  [9] tibble_2.1.3       ggplot2_3.3.0.9000 tidyverse_1.3.0    here_0.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.12        haven_2.2.0      lattice_0.20-38 
    ##  [5] snakecase_0.11.0 colorspace_1.4-1 vctrs_0.2.2      generics_0.0.2  
    ##  [9] htmltools_0.4.0  yaml_2.2.0       rlang_0.4.4      pillar_1.4.3    
    ## [13] withr_2.1.2      glue_1.3.1       DBI_1.1.0        dbplyr_1.4.2    
    ## [17] modelr_0.1.5     readxl_1.3.1     lifecycle_0.1.0  munsell_0.5.0   
    ## [21] gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5      evaluate_0.14   
    ## [25] knitr_1.27       fansi_0.4.1      broom_0.5.3      Rcpp_1.0.3      
    ## [29] backports_1.1.5  scales_1.1.0     jsonlite_1.6     fs_1.3.1        
    ## [33] hms_0.5.3        digest_0.6.23    stringi_1.4.5    grid_3.6.1      
    ## [37] rprojroot_1.3-2  cli_2.0.1        tools_3.6.1      magrittr_1.5    
    ## [41] crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0   xml2_1.2.2      
    ## [45] reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1 rmarkdown_2.1   
    ## [49] httr_1.4.1       rstudioapi_0.10  R6_2.4.1         nlme_3.1-143    
    ## [53] compiler_3.6.1
