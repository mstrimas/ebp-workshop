library(auk)
library(here)
library(tidyverse)

# countries
f_ebd <- here("raw-data", "ebd_relApr-2019_yucatan.txt")
f_sed <- here("raw-data", "ebd_sampling_relApr-2019_yucatan.txt")
auk_ebd("ebd_relApr-2019.txt", "ebd_sampling_relApr-2019.txt") %>% 
  auk_country(c("Mexico", "Guatemala", "Belize")) %>% 
  auk_date(c("2014-01-01", "2018-12-31")) %>% 
  auk_filter(f_ebd, f_sed)

# mexican states in yucatan
f_ebd_workshop <- here("raw-data", "ebd_2015-2016_yucatan.txt")
f_sed_workshop <- here("raw-data", "ebd_sampling_2015-2016_yucatan.txt")
auk_ebd(f_ebd, f_sed) %>% 
  auk_state(c("MX-CAM", "MX-CHP", "MX-TAB", "MX-ROO", "MX-YUC")) %>% 
  auk_date(c("2015-01-01", "2016-12-31")) %>% 
  auk_filter(f_ebd_workshop, f_sed_workshop)
# other countries
f_ebd_nomx <- here("raw-data", "ebd_2015-2016_nomx.txt")
f_sed_nomx <- here("raw-data", "ebd_sampling_2015-2016_nomx.txt")
auk_ebd(f_ebd, f_sed) %>% 
  auk_country(c("Guatemala", "Belize")) %>% 
  auk_date(c("2015-01-01", "2016-12-31")) %>% 
  auk_filter(f_ebd_nomx, f_sed_nomx)

# trim headers and combine
str_glue("sed '1d' {f_ebd_nomx} >> {f_ebd_workshop}") %>% system()
str_glue("sed '1d' {f_sed_nomx} >> {f_sed_workshop}") %>% system()

# create a basic download file
zip(here("raw-data", "workshop-data.zip"), c(f_ebd_workshop, f_sed_workshop))

# zip individually
str_glue("gzip < {f_ebd_workshop} > {paste0(f_ebd_workshop, '.gz')}") %>% 
  system()
str_glue("gzip < {f_sed_workshop} > {paste0(f_sed_workshop, '.gz')}") %>% 
  system()

# create packages
str_glue("mv {paste0(f_ebd_workshop, '.gz')} {here('raw-data', 'ebd-files')}") %>% 
  system()
str_glue("tar -cf {str_replace(f_ebd_workshop, 'txt$', 'tar')} ",
         "-C {here('raw-data', 'ebd-files')} .") %>% 
  system()
str_glue("rm {here('raw-data', 'ebd-files', '*.gz')}") %>% 
  system()

str_glue("mv {paste0(f_sed_workshop, '.gz')} {here('raw-data', 'ebd-files')}") %>% 
  system()
str_glue("tar -cf {str_replace(f_sed_workshop, 'txt$', 'tar')} ",
         "-C {here('raw-data', 'ebd-files')} .") %>% 
  system()
str_glue("rm {here('raw-data', 'ebd-files', '*.gz')}") %>% 
  system()

# clean up
unlink(c(f_ebd_nomx, f_sed_nomx))