mc_samp_data <- readRDS("/home/mtr/rfactory/data_test.RDS")

mc_samp <- lapply(mc_samp_data$splits, function(x)  x$in_id)

saveRDS(mc_samp, "read_mc_samp_.RDS")
