distdiff2 <-function(d1, d2){
  v <- d1-d2
  d <- density(v)
  pdf <- cbind(d$x, d$y/sum(d$y))
  colnames(pdf) <- c("x", "y")
  return(pdf)
}

d_hab_hab <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                            |ndvi_raw$cropland == 176
                                            & ndvi_raw$NDVI_count > 5],
                              5000, replace = TRUE),
                       sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                            |ndvi_raw$cropland == 176
                                            & ndvi_raw$NDVI_count > 5],
                              5000, replace = TRUE)
)

dvar_hab_hab <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                            |ndvi_raw$cropland == 176
                                            & ndvi_raw$NDVI_count > 5],
                              5000, replace = TRUE),
                       sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                            |ndvi_raw$cropland == 176
                                            & ndvi_raw$NDVI_count > 5],
                              5000, replace = TRUE)
)

d_hab_corn <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI[ndvi_raw$cropland == 1
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
)

dvar_hab_corn <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
)

dvar_hab_wheat <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                      |ndvi_raw$cropland == 176
                                                      & ndvi_raw$NDVI_count > 5],
                                   5000, replace = TRUE),
                            sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 24
                                                      & ndvi_raw$NDVI_count > 5],
                                   5000, replace = TRUE)
)

d_hab_wheat <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                                      |ndvi_raw$cropland == 176
                                                      & ndvi_raw$NDVI_count > 5],
                                   5000, replace = TRUE),
                            sample(ndvi_raw$NDVI[ndvi_raw$cropland == 24
                                                      & ndvi_raw$NDVI_count > 5],
                                   5000, replace = TRUE)
)

dvar_hab_gum <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                    |ndvi_raw$cropland == 176
                                                    & ndvi_raw$NDVI_count > 5],
                                 5000, replace = TRUE),
                          sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 4
                                                    & ndvi_raw$NDVI_count > 5],
                                 5000, replace = TRUE)
)

d_hab_gum <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                                    |ndvi_raw$cropland == 176
                                                    & ndvi_raw$NDVI_count > 5],
                                 5000, replace = TRUE),
                          sample(ndvi_raw$NDVI[ndvi_raw$cropland == 4
                                                    & ndvi_raw$NDVI_count > 5],
                                 5000, replace = TRUE)
)

dvar_hab_bare <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                             |ndvi_raw$cropland == 176
                                             & ndvi_raw$NDVI_count > 5],
                               5000, replace = TRUE),
                        sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 61
                                             & ndvi_raw$NDVI_count > 5],
                               5000, replace = TRUE)
)

d_hab_bare <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                             |ndvi_raw$cropland == 176
                                             & ndvi_raw$NDVI_count > 5],
                               5000, replace = TRUE),
                        sample(ndvi_raw$NDVI[ndvi_raw$cropland == 61
                                             & ndvi_raw$NDVI_count > 5],
                               5000, replace = TRUE)
)

dvar_hab_alph <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 36
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
)

d_hab_alph <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI[ndvi_raw$cropland == 36
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
)

dvar_crp_crp <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                                    |ndvi_raw$cropland == 24
                                                    |ndvi_raw$cropland == 36
                                                    |ndvi_raw$cropland == 4
                                                    & ndvi_raw$NDVI_count > 5],
                                 20000, replace = TRUE),
                          sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                                    |ndvi_raw$cropland == 24
                                                    |ndvi_raw$cropland == 36
                                                    |ndvi_raw$cropland == 4
                                                    & ndvi_raw$NDVI_count > 5],
                                 20000, replace = TRUE)
)

d_crp_crp <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 1
                                                    |ndvi_raw$cropland == 24
                                                    |ndvi_raw$cropland == 36
                                                    |ndvi_raw$cropland == 4
                                                    & ndvi_raw$NDVI_count > 5],
                                 20000, replace = TRUE),
                          sample(ndvi_raw$NDVI[ndvi_raw$cropland == 1
                                                    |ndvi_raw$cropland == 24
                                                    |ndvi_raw$cropland == 36
                                                    |ndvi_raw$cropland == 4
                                                    & ndvi_raw$NDVI_count > 5],
                                 20000, replace = TRUE)
)
