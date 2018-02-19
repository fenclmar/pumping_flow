# ===========================
# flow generator functions:
# ===========================

# Functions for modeling realistic I/I patterns



# ---------------------------------
add_II_to_SPGflow <- function (spg_flow, fast = NULL, slow = NULL, gw = NULL) {
    
    # add inflow/infiltration component to the SPG flow object
    # for more infro on sewage pattern generator (SPG) see
    # https://github.com/scheidan/SPG
    #
    # Args: 
    # spg_flow - SPG flow object with sewage flow and substance concentration
    # fast - time series (xts) with fast inflow
    # slow - time series (xts) with slow inflow
    # slow - time series (xts) with groundwater infiltration
    
    require(xts)
    
    
    Q <- merge(slow, fast)
    Q <- merge(Q, gw)
    
    
    st_sim <- as.POSIXct(as.Date(index(Q)[1]))  # beggining of simulation to 00:00:00
    
    
    w_index <- seq(st_sim, length = nrow(flow_s1),
                   by = attributes(flow_s1)$temp.res.sim)
    Q_waste <- xts(flow_s1[ ,1], w_index)
    S_waste <- xts(flow_s1[ ,2], w_index)
    Q <- merge(Q, Q_waste)
    Q <- merge(Q, S_waste)
    
    
    # treat NA values (different time steps)
    for (i in 1 : ncol(Q)) {
        Q[, i] <- na.locf(Q[, i])    
    }
    
    if( is.regular(Q) ) {
        new_resolution <- difftime(index(Q[2]), index(Q[1]), units= 'secs')
    } else {
        new_resolution <- NA
        warning('resulting time series is not regular!')
    }
    
    
    Q_tot <- apply(Q[ , 1:(ncol(Q) - 1)], 1, sum, na.rm = T)
    flow_new <- as.matrix(data.frame('Q' = Q_tot, 'S'= Q$S_waste))
    
    # Add I/I runoff to wastewater flow object
    att <- attributes(spg_flow)
    spg_flow <- flow_new
    attributes(spg_flow)$class <- 'flow'
    attributes(spg_flow)$temp.res.sim <- as.numeric(new_resolution)
    
    wd_we <- get_weekday_or_weekend(index(Q))
    wd_we <- paste(1 : length(wd_we), wd_we, sep = "_")
    
    attributes(spg_flow)$dimnames[[1]] <- wd_we
    
    attributes(spg_flow)$timeindex <- index(Q)
    attributes(spg_flow)$wastewater <- spg_flow[ ,1] 
    if ( !is.null(fast) ) { attributes(spg_flow)$fastrunoff <- Q$fast }
    if ( !is.null(slow) ) { attributes(spg_flow)$slowrunoff <- Q$slow }
    if ( !is.null(gw) ) { attributes(spg_flow)$groundwater <- Q$gw }
    
    
    return(spg_flow)
}



# ---------------------------------
get_weekday_or_weekend <- function (tim) {
    
    # arg: tim - vector with POSIXtimes
    # output: wd_we - vector corresponding to POSIXtimes indicating if
    #                 the time is weekday ('WD') or Weekend ('WE')
    
    day_names <- weekdays(tim, abbreviate = T)
    day_ids <- match(day_names, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
    
    wd_we <- rep('WD', length(day_ids))
    wd_we[ day_ids > 5 ] <- 'WE'
    
    return(wd_we)
}
# ---------------------------------
