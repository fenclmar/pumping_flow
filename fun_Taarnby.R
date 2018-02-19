# ===========================
# level to volume functions
# ===========================

# - Functions to convert level to volume at Taarnby catchment

# ---------------------------------
level_to_volume_p8_10 <- function (h) {
    
    # calculate volume from level for pump sumps PV_8 - 10.
    # h - pump level
    
    ha_tab0 <- matrix(c(0.0, 0.3701522192,
                        0.410, 2.679839769,
                        0.684, 4.578881589,
                        1.139, 9.701829835,
                        1.641, 12.22696363,
                        5.500, 12.22696363), 6, 2, byrow = T)
    
    
    if(h < 0) { stop ('level has to be positive number!')}
    
    # interpolate to get area at level h
    a_h <- approx(ha_tab0[ ,1], ha_tab0[ ,2], xout = h)$y
    ha_tab <- rbind(ha_tab0, c(h, a_h))
    ha_tab <- ha_tab[order(ha_tab[ ,1]), ]
    ha_tab <- ha_tab[which(ha_tab[ ,1] <= h), ]
    
    # calculate volume
    s <- ha_tab[-1, 1] - ha_tab[-nrow(ha_tab), 1]
    h_mid <- ha_tab[-1, 1] - s/2
    a_mid <- approx(ha_tab0[ ,1], ha_tab0[ ,2], xout = h_mid)$y
    vol <- sum(s * a_mid)
    
    return(vol)
    
}
# ---------------------------------



# ---------------------------------
level_to_volume_pi1 <-  function (h) {
    
    # calculate volume from level for pump sump PVI_1.
    # h - pump level
    
    if (length(h[h < 0]) > 0) { stop ('level has to be positive number!')}
    if (length(h[h > 4.02]) > 0) {
        warning (paste(length(h[h > 4.02]),
                       'values of water level exceeded pump sump height!'))
    } 
    
    return(pi * 2.8 ^ 2 * h / 4)
}
# ---------------------------------



# ---------------------------------
level_to_volume <- function (h, p_name){
    
    # get volume
    
    p_id <- match(p_name, c('PVI_1', 'PV_8', 'PV_9', 'PV_10'))
    
    
    if (is.na(p_id) == T){
        stop('pump name does not match Tarnby pumps')
    } else {
        if(p_id == 1) {
            return(level_to_volume_pi1(h))
        } else {
            return(level_to_volume_p8_10(h))
        }
        
    }
    
}
# ---------------------------------

