library(tidyverse)
library(here)
library(lubridate)

mg <- read_tsv(here::here("consolidatR", "tara_simons_timeseries_basic_meta.tsv")) 

### formatting time series
timeseries.isdcm <- mg %>% 
  filter(sample_category=="timeseries") %>%
  mutate(environment=case_when(depth < 150 & depth > 25 ~ "dcm",
                               depth <= 25 ~ "surface",
                               depth > 150 & depth < 300 ~ "mixed.layer",
                               depth > 300 ~ "mesopelagic")) %>%
  mutate(is.dcm=ifelse(environment=="dcm", "yes", "no")) %>%
  select(-depth_category, -replicate, -time, -station_num_BODC, -bottle_num_BODC) %>%
  mutate(date = ymd(date))

write_tsv(timeseries.isdcm, here::here("results", "timeseries_metadata.tsv"))
write_rds(timeseries.isdcm, here::here("results", "timeseries_metadata"))

### formatting TARA oceans data
tara.1 <- read_tsv(here("consolidatR", "tara_metadata_slh_v1.tsv"))
  
tara.0 <- read_tsv(here("consolidatR", "tara_env_feat.tsv"))  %>%
  separate(run_acc, into=c("a", "b", "c", "d", "e", "f", "g"), remove = T, sep = "\\|") %>%
  pivot_longer(cols = c("a", "b", "c", "d", "e", "f", "g"), names_to = "name", values_to = "ebi_run_id") %>%
  rename(ebi_sample_id=sample_acc) %>%
  select(-name) %>%
  filter(!is.na(ebi_run_id))
  
tara.final <- left_join(tara.1, tara.0) %>%
  mutate(a=str_extract(env_feat, "^\\(\\w+\\)")) %>%
  mutate(new_env=case_when(str_detect(a, "DCM") ~ "dcm",
                           str_detect(a, "SRF") ~ "surface",
                           str_detect(a, "MIX") ~ "mixed.layer",
                           str_detect(a, "MES") ~ "mesopelagic",
                           is.na(a) ~ environment)) %>%
  mutate(new_env=ifelse(new_env == "srf", "surface", new_env)) %>%
  select(-environment, -a, -env_feat, -sequencing_type, -sample_id) %>%
  rename(environment=new_env, sampleID=lib_id, lat=latitude, lon=longitude) %>%
  mutate(is.dcm=ifelse(environment=="dcm", "yes", "no")) %>%
  mutate(date = mdy(date),
         sample_category="tara")

write_tsv(tara.final, here::here("results", "tara_metadata.tsv"))
write_rds(tara.final, here::here("results", "tara_metadata"))

## GEOTRACES DCM
GA02.dcm <- read_tsv(here::here("CTD2mg", "GA02_dcm.tsv")) %>% mutate(section="GA02")
GA03.dcm <- read_tsv(here::here("CTD2mg", "GA03_dcm.tsv")) %>% mutate(section="GA03")
GA10.dcm <- read_tsv(here::here("CTD2mg", "GA10_dcm.tsv")) %>% mutate(section="GA10")
GP13.dcm <- read_tsv(here::here("CTD2mg", "GP13_dcm.tsv")) %>% mutate(section="GP13")

dcm <- bind_rows(GA02.dcm, GA03.dcm, GA10.dcm, GP13.dcm) %>%
  group_by(sampleID) %>%
  arrange(dcm.depth.diff) %>%
  slice(1)

write_tsv(dcm, here::here("results", "geotraces_dcm.tsv"))
write_rds(dcm, here::here("results", "geotraces_dcm"))

## GEOTRACES biogeochemistry

GA02.btl <- read_tsv(here::here("bottle2mg", "GA02_mg_bottles_long.tsv")) %>% mutate(section="GA02")
GA03.btl <- read_tsv(here::here("bottle2mg", "GA03_mg_bottles_long.tsv")) %>% mutate(section="GA03")
GA10.btl <- read_tsv(here::here("bottle2mg", "GA10_mg_bottles_long.tsv")) %>% mutate(section="GA10")
GP13.btl <- read_tsv(here::here("bottle2mg", "GP13_mg_bottles_long.tsv")) %>% mutate(section="GP13")

geo <- bind_rows(GA02.btl, GA03.btl, GA10.btl, GP13.btl)

write_tsv(geo, here::here("results", "geotraces_biogeochem.tsv"))
write_rds(geo, here::here("results", "geotraces_biogeochem"))
