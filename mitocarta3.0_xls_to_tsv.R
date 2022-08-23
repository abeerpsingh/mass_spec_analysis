library(readxl)
library(readr)
mitocarta <- read_xls(path = "Human.MitoCarta3.0.xls",
                      sheet = 2) %>%
  as_tibble()
write_tsv(x = mitocarta, file = "mitocarta3.0.tsv")
