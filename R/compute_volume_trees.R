compute_volume_trees <- function(data_trees) {
  
  data_trees %>% 
    
    dplyr::mutate(
      # Volumes of the upper fourths
      c_top = h_m - hmax_m,
      v_upper_n = (4/3) * pi * rn_m^2 * c_top / 8,
      v_upper_s = (4/3) * pi * rs_m^2 * c_top / 8,
      v_upper_w = (4/3) * pi * rw_m^2 * c_top / 8,
      v_upper_e = (4/3) * pi * re_m^2 * c_top / 8,
      v_upper = v_upper_n + v_upper_s + v_upper_w + v_upper_e,
      
      # Volumes of the lower fourths
      c_bot = hmax_m - hbase_m,
      v_lower_n = (4/3) * pi * rn_m^2 * c_bot / 8,
      v_lower_s = (4/3) * pi * rs_m^2 * c_bot / 8,
      v_lower_w = (4/3) * pi * rw_m^2 * c_bot / 8,
      v_lower_e = (4/3) * pi * re_m^2 * c_bot / 8,
      v_lower = v_lower_n + v_lower_s + v_lower_w + v_lower_e,
      
      # Volume of the assymetric crown
      volume_m3 = v_upper + v_lower
    )
  
}
