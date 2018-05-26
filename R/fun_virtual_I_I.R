# ===========================
# flow generator functions:
# ===========================

# Functions for modeling realistic I/I patterns


# ---------------------------------
simulate_fast_runoff <- function(rain, area_fast, reserv_count = 3,
                                 retention_time, temp_res) {
    # Arguments:
    # rain - rainfall intensity (mm/h)
    # area_fast - fast runoff area (ha)
    # reserv_count - reservoir count
    # retention time - retention time (s)
    # temp_res - temporal resolution of generated flow time series (s)
    # Outputs:
    # Q_fast - fast runoff xts series (l/s)
    
    require(Rcpp)
    sourceCpp("muskingum_simulation.cpp")
    #sourceCpp("MuskingumReservoirsInSeries.cpp")
    
    flow_index <- seq(from = start(rain), to = end(rain), by = temp_res)
    xts_flow <- zoo(flow_index, flow_index)
    
    
    flow_fast <- muskingum_reservoirs_simulation(reserv_count, retention_time,
                                                 rain, xts_flow)
    flow_fast <- xts(flow_fast, flow_index)
    Q_fast <- flow_fast * area_fast * 10 / 3.6 # Q fast [l/s]
    
    return(Q_fast)
    
}



# ---------------------------------
simulate_slow_runoff <- function(rain, area_slow, reserv_count = 3, 
                                 retention_time, temp_res) {
                 
    # Arguments:
    # rain - rainfall intensity (mm/h)
    # area_slow - slow runoff area (ha)
    # reserv_count - reservoir count
    # retention time - retention time (s)
    # temp_res - temporal resolution of generated flow time series (s)
    # Outputs:
    # Q_slow - slow runoff xts series (l/s)
    
    require(Rcpp)
    sourceCpp("muskingum_simulation.cpp")
    #sourceCpp("MuskingumReservoirsInSeries.cpp")
    
    flow_index <- seq(from = start(rain), to = end(rain), by = temp_res)
    xts_flow <- zoo(flow_index, flow_index)
    
    flow_slow <- muskingum_reservoirs_simulation(reserv_count, retention_time,
                                                 rain, xts_flow)
    flow_slow <- xts(flow_slow, flow_index)
    Q_slow <- flow_slow * area_slow * 10 / 3.6 # Q slow [l/s]

    return(Q_slow)
    
}


# ---------------------------------
add_II <- function (rain, flow_sewerage, area_fast, area_slow,
                         reserv_count_fast = 3, reserv_count_slow = 3,
                         retention_time_fast = 9000, retention_time_slow = 60000,
                         flow_gw = NULL, temp_res) {
    
    # generate I_I  and create an object compatible with SPG.
    # The object containins flow components and rainfall (used to generate I_I)

    flow_fast <- simulate_fast_runoff(rain, area_fast, reserv_count_fast,
                                      retention_time_fast, temp_res)
    
    flow_slow <- simulate_slow_runoff(rain, area_slow, reserv_count_slow,
                                      retention_time_slow, temp_res)
    
    total_runoff <- add_II_to_SPGflow (flow_sewerage, flow_fast, flow_slow,
                                       flow_gw)

    # add rainfall info
    
    attributes(total_runoff)$rain <- rain
    
    return(total_runoff)
}


# ---------------------------------
add_II_to_SPGflow <- function (spg_flow, fast = NULL, slow = NULL, gw = NULL) {
    
    # add inflow/infiltration component to the SPG flow object
    # for more info on sewage pattern generator (SPG) see
    # https://github.com/scheidan/SPG
    #
    # Args: 
    # spg_flow - SPG flow object with sewage flow and substance concentration
    # fast - time series (xts) with fast inflow
    # slow - time series (xts) with slow inflow
    # gw - time series (xts) with groundwater infiltration
    
    require(xts)
    
    Q <- merge(slow, fast)
    Q <- merge(Q, gw)

    st_sim <- as.POSIXct(as.Date(index(Q)[1])) + 60 # beggining of simulation to 00:00:01
    
    
    w_index <- seq(st_sim, length = nrow(spg_flow),
                   by = attributes(spg_flow)$temp.res.sim)
    
    Q_waste <- xts(spg_flow[ ,1], w_index)
    S_waste <- xts(spg_flow[ ,2], w_index)
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
simulate_pumping <- function(pump, inflow_subc, pumping_up = NULL) {
    # generate virtual pumping possibly including as an attribute information on pumping from from a singled upstream pump
    # Arguments:
    # pump - SPG pump object
    # inflow_subc - inflow from the catchment (spg flow object)
    # pumping_up - inflow from the upstream pump (spg flow object)

    if (is.null(pumping_up)) {
        ps_flow <- pump (inflow_subc)
        attributes(ps_flow)$inflow <- inflow_subc[ ,1]
    } else {
        ps_flow <- pump (inflow_subc + pumping_up)
        attributes(ps_flow)$inflow <- (inflow_subc + pumping_up)[ ,1]
        attributes(ps_flow)$pump.state.up <- attributes(pumping_up)$pump.state
        attributes(ps_flow)$inflow.up <- pumping_up[ ,1]
    }
    

    return(ps_flow)
}



# ---------------------------------
convert_pumping_to_xts <- function (ps_flow, time_index,
                                    inflow = T, upstream_inflow = T) {
    
    # generate xts time series with multiple records of pumping station
    # Arguments:
    # ps_flow - flow object (supplied e.g. by gsimulate_pumping function)
    # time_index - vector with POSIX time. It has to match number of ps_flow
    #              records
    # Outputs:
    # pump_data - xts time series with pumping station data
    #             (in a similar from as usually provided by SCADA system)

    require(xts)
    
    if (length(time_index) != nrow(ps_flow)) {
        stop ('Number of ps rows do not match time_index length!')
    }
    require(xts)
    
    pump_data <- data.frame('time' = time_index)
    pump_data$vol <- attributes(ps_flow)$V.sump
    pump_data$cur <- as.numeric(attr(ps_flow, 'pump.state') == 'on')
    
    if (!is.null(attributes(ps_flow)$pump.state.up)) {
        pump_data$curup <- attr(ps_flow, 'pump.state.up') == 'on'
    }
    
    
    if (inflow == T){
        pump_data$Qin <- attr(ps_flow, 'inflow')
    }
    

    if (upstream_inflow == T) {
        if (!is.null(attributes(ps_flow)$pump.state.up)) {    
            pump_data$Qup <- attr(ps_flow, 'inflow.up')
        } else {
            warning('There is no upstream inflow!')
        }
    }

    
    pump_data <- xts(pump_data[ ,-1], pump_data[, 1])

    return(pump_data)

}



# ---------------------------------
sample_pump_data <- function (pump_data, freq, cur_discreate = T) {
    
    # Sample from xts pumping station time series to simulate raw sampling
    # Arguments:
    # pump_data - xts time series with pumping station data
    #             (in a similar form as usually provided by SCADA system)
    # freq - frequency o sampling
    # cur_continous - does current record represent discreate value at one given 
    #                 time (T), or it is proportional to the duration of pumping
    #                 within the time step
    #                 
    # Outputs:
    # pump_resampled - xts time series with resampled pumping station data
    require(xts)
        
    # sample volume
    vol <- pump_data$vol
    vol_resampled <- vol[round(seq(1/freq, nrow(pump_data), 1/freq),1)]    
    
    # sample current
    cur_columns <- which(is.element(colnames(pump_data), c('cur', 'curup')))
    cur <- pump_data[ , cur_columns]

    pump_resampled <- pump_data[round(seq(1/freq, nrow(pump_data), 1/freq),1), ]    
    

    if (cur_discreate) {
        cur_resampled <- cur[round(seq(1/freq, nrow(pump_data), 1/freq),1), ]
    } else {
        cur_resampled <- rollapply(as.zoo(cur), width = 1/freq, mean,
                                   by = 1/freq, align = 'right')
    }
    
    # sample inflow
    Q_columns <- which(is.element(colnames(pump_data), c('Qin', 'Qup')))
    Q <- pump_data[ , Q_columns]
    Q_resampled <- rollapply(as.zoo(Q), width = 1/freq, mean,
                             by = 1/freq, align = 'right')
    
    pump_resampled <- merge(vol_resampled, cur_resampled, Q_resampled)
    colnames(pump_resampled) <- colnames(pump_data)

    
    return(pump_resampled)
}

