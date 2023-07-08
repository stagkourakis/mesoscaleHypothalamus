# Plot probes in Mouse Common Coordinate Framework (CCF) v3 space
#
# From Allen Brain: If you use cocoframer to make your own 3D brain structure
# images and animations, please cite:
#   © 2018 Allen Institute for Brain Science. Allen Brain Explorer. Available
#   from: connectivity.brain-map.org/3d-viewer/

library(cocoframer)
library(fs)
library(htmlwidgets)
library(magrittr)
library(rgl)
library(tidyverse)

# At this point you need to specify the directory where you have placed your
# probe coordinates;
# go to Session - Set Working Directory - etc.
# C:\Users\Bumblebee\Desktop\Ranalysis


# probe coordinates -------------------------------------------------------

files <- dir_ls() %>% keep(~path_ext(.x) == "csv")

ids <- files %>% str_extract("\\d+")

coords <- 
  files %>%
  map(read_csv) %>%
  map(rename_all, str_to_lower) %>%
  map2_dfr(ids, ~mutate(.x, id = .y))

# plot functions ----------------------------------------------------------

plot_probe_pts <- function(coords, color) {
  
  # y and z axes are inverted in CCF v3
  # compare plot of ref brain using axes3d() vs. coords in Fiji
  
  coords_invert <- coords %>% mutate(y = 8000 - y, z = 11375 - z)
  
  pch3d(
    x = coords_invert$x,
    y = coords_invert$y,
    z = coords_invert$z,
    pch = 20,
    cex = 1,
    radius = 300,
    color = color
  )
}

plot_probe_line <- function(coords, color) {
  
  # ref: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
  
  coord_means <-
    coords %>%
    select(x, y, z) %>%
    summarise(across(everything(), mean))

  coords_centered <-
    coords %>%
    transmute(
      x = x - coord_means$x,
      y = y - coord_means$y,
      z = z - coord_means$z
    )
  
  pca_comp_1 <-
    coords_centered %>%
    as.matrix() %>%
    princomp() %$%
    loadings %>%
    magrittr::extract(, 1)

  # sequence of points along first component; ensure sequence encompasses
  # probe's start/end points and step size is appropriate
  
  probe_pts_range <-
    seq(-30000, 30000, 10) %>%
    as.matrix() %>%
    t()
  
  probe_pts_line <-
    t(pca_comp_1 %*% probe_pts_range) %>%
    set_colnames(c("x", "y", "z")) %>%
    as_tibble() %>%
    mutate(
      x = x + coord_means$x,
      y = y + coord_means$y,
      z = z + coord_means$z
    )

 
  # get point along first component at the probe's most ventral point (i.e. the
  # max point along the y axis); assumes the most ventral coordinate recorded is
  # the most ventral point of the probe
  
  coords_y_max <- coords %>% pull(y) %>% max()
  
  probe_pt_ventral <-
    probe_pts_line %>%
    mutate(y_dist_from_max = abs(y - coords_y_max)) %>%
    filter(y_dist_from_max == min(y_dist_from_max))

  # get point along first components at approx. most dorsal point of the
  # reference brain
  
  probe_pt_dorsal <-
    probe_pts_line %>%
    mutate(y_dist_from_min = abs(y - 300)) %>%
    filter(y_dist_from_min == min(y_dist_from_min))
  
  # y and z axes are inverted in CCF v3
  # compare plot of ref brain using axes3d() vs. coords in Fiji
  
  probe_line <-
    bind_rows(probe_pt_dorsal, probe_pt_ventral) %>%
    select(x, y, z) %>%
    mutate(y = 8000 - y, z = 11375 - z)

  lines3d(probe_line, color = color, lwd = 4.0)
}

# plots -------------------------------------------------------------------

structures <- c("root", "HY")

mesh_list <- map(structures, ccf_2017_mesh)

names(mesh_list) <- structures

open3d()

plot_ccf_meshes(
  mesh_list,
  fg_structure = "HY", 
  fg_alpha = 0.4,
  bg_structure = "root",
  bg_alpha = 0.4
)

# ref: https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=4
id_colors <-
  # c('#1b9e77', '#d95f02', '#7570b3', '#e7298a') %>%
  c('#ff00ff', '#000000', '#ff0000','#0000ff', '#00ff00', '#754096') %>%
  map(rep, 8) %>%
  flatten_chr()

# plot probes
coords %>%
  split(list(coords$probe, coords$id)) %>%
  walk2(id_colors, ~plot_probe_line(.x, color = .y))

# plot coordinates
coords %>%
  split(list(coords$probe, coords$id)) %>%
  walk2(id_colors, ~plot_probe_pts(.x, color = .y))

# save 3D image
plot_widget <- rglwidget(
  scene3d(),
  width = 600, 
  height = 600
)

saveWidget(plot_widget, "2021-01-30_probes.html")

# ─ Session info ────────────────────────────────────────────────────────────
# setting  value                       
# version  R version 3.5.1 (2018-07-02)
# os       macOS  10.15.7              
# system   x86_64, darwin15.6.0        
# ui       RStudio                     
# language (EN)                        
# collate  en_US.UTF-8                 
# ctype    en_US.UTF-8                 
# tz       America/New_York            
# date     2021-01-30                  
# 
# ─ Packages ────────────────────────────────────────────────────────────────
# package          * version  date       lib source                                    
# assertthat         0.2.1    2019-03-21 [1] CRAN (R 3.5.1)                            
# backports          1.1.5    2019-10-02 [1] CRAN (R 3.5.2)                            
# broom              0.5.5    2020-02-29 [1] CRAN (R 3.5.2)                            
# callr              3.5.1    2020-10-13 [1] CRAN (R 3.5.1)                            
# cellranger         1.1.0    2016-07-27 [1] CRAN (R 3.5.0)                            
# cli                2.2.0    2020-11-20 [1] CRAN (R 3.5.1)                            
# cocoframer       * 0.1.1    2020-12-15 [1] Github (AllenInstitute/cocoframer@1de30a8)
# colorspace         2.0-0    2020-11-11 [1] CRAN (R 3.5.1)                            
# crayon             1.3.4    2017-09-16 [1] CRAN (R 3.5.0)                            
# crosstalk          1.0.0    2016-12-21 [1] CRAN (R 3.5.0)                            
# DBI                1.1.0    2019-12-15 [1] CRAN (R 3.5.2)                            
# dbplyr             1.4.2    2019-06-17 [1] CRAN (R 3.5.2)                            
# desc               1.2.0    2018-05-01 [1] CRAN (R 3.5.0)                            
# devtools           2.2.2    2020-02-17 [1] CRAN (R 3.5.2)                            
# digest             0.6.27   2020-10-24 [1] CRAN (R 3.5.1)                            
# dplyr            * 1.0.2    2020-08-18 [1] CRAN (R 3.5.1)                            
# ellipsis           0.3.1    2020-05-15 [1] CRAN (R 3.5.1)                            
# evaluate           0.14     2019-05-28 [1] CRAN (R 3.5.2)                            
# fansi              0.4.1    2020-01-08 [1] CRAN (R 3.5.2)                            
# fastmap            1.0.1    2019-10-08 [1] CRAN (R 3.5.2)                            
# forcats          * 0.5.0    2020-03-01 [1] CRAN (R 3.5.2)                            
# fs               * 1.3.2    2020-03-05 [1] CRAN (R 3.5.2)                            
# generics           0.1.0    2020-10-31 [1] CRAN (R 3.5.1)                            
# ggplot2          * 3.3.2    2020-06-19 [1] CRAN (R 3.5.1)                            
# glue               1.4.2    2020-08-27 [1] CRAN (R 3.5.1)                            
# gtable             0.3.0    2019-03-25 [1] CRAN (R 3.5.1)                            
# haven              2.2.0    2019-11-08 [1] CRAN (R 3.5.2)                            
# hms                0.5.3    2020-01-08 [1] CRAN (R 3.5.2)                            
# htmltools          0.5.0    2020-06-16 [1] CRAN (R 3.5.1)                            
# htmlwidgets      * 1.5.3    2020-12-10 [1] CRAN (R 3.5.1)                            
# httpuv             1.5.2    2019-09-11 [1] CRAN (R 3.5.2)                            
# httr               1.4.2    2020-07-20 [1] CRAN (R 3.5.1)                            
# jsonlite           1.7.2    2020-12-09 [1] CRAN (R 3.5.1)                            
# knitr              1.30     2020-09-22 [1] CRAN (R 3.5.1)                            
# later              1.0.0    2019-10-04 [1] CRAN (R 3.5.2)                            
# lattice            0.20-40  2020-02-19 [1] CRAN (R 3.5.2)                            
# lifecycle          0.2.0    2020-03-06 [1] CRAN (R 3.5.2)                            
# lubridate          1.7.4    2018-04-11 [1] CRAN (R 3.5.0)                            
# magrittr         * 2.0.1    2020-11-17 [1] CRAN (R 3.5.1)                            
# manipulateWidget   0.10.1   2020-02-24 [1] CRAN (R 3.5.2)                            
# memoise            1.1.0    2017-04-21 [1] CRAN (R 3.5.0)                            
# mime               0.9      2020-02-04 [1] CRAN (R 3.5.2)                            
# miniUI             0.1.1.1  2018-05-18 [1] CRAN (R 3.5.0)                            
# modelr             0.1.6    2020-02-22 [1] CRAN (R 3.5.2)                            
# munsell            0.5.0    2018-06-12 [1] CRAN (R 3.5.0)                            
# nlme               3.1-145  2020-03-04 [1] CRAN (R 3.5.2)                            
# pillar             1.4.7    2020-11-20 [1] CRAN (R 3.5.1)                            
# pkgbuild           1.1.0    2020-07-13 [1] CRAN (R 3.5.1)                            
# pkgconfig          2.0.3    2019-09-22 [1] CRAN (R 3.5.2)                            
# pkgload            1.1.0    2020-05-29 [1] CRAN (R 3.5.1)                            
# prettyunits        1.1.1    2020-01-24 [1] CRAN (R 3.5.2)                            
# processx           3.4.5    2020-11-30 [1] CRAN (R 3.5.1)                            
# promises           1.1.0    2019-10-04 [1] CRAN (R 3.5.2)                            
# ps                 1.5.0    2020-12-05 [1] CRAN (R 3.5.1)                            
# purrr            * 0.3.4    2020-04-17 [1] CRAN (R 3.5.1)                            
# R6                 2.5.0    2020-10-28 [1] CRAN (R 3.5.1)                            
# Rcpp               1.0.5    2020-07-06 [1] CRAN (R 3.5.1)                            
# readr            * 1.3.1    2018-12-21 [1] CRAN (R 3.5.0)                            
# readxl             1.3.1    2019-03-13 [1] CRAN (R 3.5.2)                            
# remotes            2.1.1    2020-02-15 [1] CRAN (R 3.5.2)                            
# reprex             0.3.0    2019-05-16 [1] CRAN (R 3.5.2)                            
# rgl              * 0.100.50 2020-02-21 [1] CRAN (R 3.5.2)                            
# rlang              0.4.9    2020-11-26 [1] CRAN (R 3.5.1)                            
# rmarkdown          2.6      2020-12-14 [1] CRAN (R 3.5.1)                            
# rprojroot          2.0.2    2020-11-15 [1] CRAN (R 3.5.1)                            
# rstudioapi         0.13     2020-11-12 [1] CRAN (R 3.5.1)                            
# rvest              0.3.5    2019-11-08 [1] CRAN (R 3.5.2)                            
# scales             1.1.1    2020-05-11 [1] CRAN (R 3.5.1)                            
# sessioninfo        1.1.1    2018-11-05 [1] CRAN (R 3.5.0)                            
# shiny              1.4.0    2019-10-10 [1] CRAN (R 3.5.2)                            
# stringi            1.5.3    2020-09-09 [1] CRAN (R 3.5.1)                            
# stringr          * 1.4.0    2019-02-10 [1] CRAN (R 3.5.2)                            
# testthat           3.0.0    2020-10-31 [1] CRAN (R 3.5.1)                            
# tibble           * 3.0.4    2020-10-12 [1] CRAN (R 3.5.1)                            
# tidyr            * 1.0.2    2020-01-24 [1] CRAN (R 3.5.2)                            
# tidyselect         1.1.0    2020-05-11 [1] CRAN (R 3.5.1)                            
# tidyverse        * 1.3.0    2019-11-21 [1] CRAN (R 3.5.2)                            
# usethis            1.5.1    2019-07-04 [1] CRAN (R 3.5.2)                            
# utf8               1.1.4    2018-05-24 [1] CRAN (R 3.5.0)                            
# vctrs              0.3.5    2020-11-17 [1] CRAN (R 3.5.1)                            
# webshot            0.5.2    2019-11-22 [1] CRAN (R 3.5.2)                            
# withr              2.3.0    2020-09-22 [1] CRAN (R 3.5.1)                            
# xfun               0.19     2020-10-30 [1] CRAN (R 3.5.1)                            
# xml2               1.2.2    2019-08-09 [1] CRAN (R 3.5.2)                            
# xtable             1.8-4    2019-04-21 [1] CRAN (R 3.5.2)                            
# yaml               2.2.1    2020-02-01 [1] CRAN (R 3.5.2)                            
# 
# [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library