create_samsalight_inventories <- function(trees_list, exp_design) {
  
  inv_ids <- unique(exp_design$inv_id)
  invs_list <- setNames(vector("list", length(inv_ids)), inv_ids)
  
  for (id in inv_ids) {

    simu_info <- exp_design %>% 
      dplyr::filter(inv_id == id)
    
    inv_name <- unique(simu_info$inv_name)
    
    invs_list[[id]] <- trees_list[[inv_name]] %>%
      dplyr::mutate(rn_m = r_m,
                    re_m = r_m,
                    rs_m = r_m,
                    rw_m = r_m,
                    crown_type = case_match(species,
                                            c("Picea abies", "Abies alba") ~ "P",
                                            "Fagus sylvatica" ~ "E"),
                    hmax_m = case_match(crown_type,
                                        "P" ~ hbase_m,
                                        "E" ~ hbase_m + 1/2*(h_m - hbase_m))) %>%
      dplyr::select(-r_m, -crown_lad) %>% 
      dplyr::left_join(simu_info %>% dplyr::select(sp, crown_lad = lad),
                       by = c("species" = "sp"))
  }
  
  invs_list
}