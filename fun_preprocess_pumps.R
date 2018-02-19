# ===============================
# Preprocessing of pumping data
# ===============================

# Preprocess pumping data, i.e. treat irregular time step, missing records and
# outliers


# ---------------------------------
round_zoo_to <- function (z, n_seconds) {
    
    # arg: z - zoo otime series
    #      n_seconds - time in seconds to round time stamps to 
    # out: zoo with rounded indexes
    
    require(zoo)
    
    t <- as.numeric(index(z))
    t <- round(t / n_seconds, 0) * n_seconds
    t <- as.POSIXct(t, origin = '1970-01-01 00:00:00')
    
    z2 <- zoo(coredata(z), t)
    
    return(z2)
}
# ---------------------------------



# ---------------------------------
filter_redundant_obs <- function (z, t_tsh) {
    
    # Function for filtering redundant data points
    # Parameters: z - zoo time series (vector)
    #             t_ths - max time interval between two points [seconds] to
    #                     asssume second of them as redundant
    
    require(zoo)
    
    t_s <- index(z)
    t_s <- as.numeric(z)
    step <- t_s[-1] - t_s[-length(t_s)]
    id_a <- which (step < t_tsh)  # indexes of values where time steps to the
    # next observation is shorter than tsh
    id_b <- id_a + 1 
    
    while (length(id_a) > 0){
        #print(length(id_a))
        id_b <- id_a + 1  # index of observation closer than 100 s to the
        # previous one 
        z <- z[-id_b]
        
        step <- as.numeric(index(z)[-1]) - as.numeric(index(z)[-length(z)])
        id_a <- which (step < t_tsh)  # indexes of values where time s
    }
    
    return(z)
}   
# ---------------------------------



# ---------------------------------
insert_times <- function (st, en, mean_int) {
    
    # generate time stamps between two times with regular spacing that 
    # spacing is close to mean_interval
    # Arguments: st - POSIX time
    #            en - POSIX time
    #            mean_int - mean interval
    
    gap <- as.numeric(en) - as.numeric(st)
    
    if (gap < mean_int * 4 / 3) {
        return (st + numeric(0))
    } else {
        n <- floor(gap / mean_int)
        e1 <- abs(mean_int - gap / (n + 1))
        e2 <- abs(mean_int - gap / (n + 2))
        
        x1 <- seq(st, en, length = n + 1)
        x1 <- seq(st, en, length = n + 2)
        
        if (e1 <= e2) {
            return(x1[- c(1, length(x1))])
        } else {
            return(x2[- c(1, length(x2))])
        }
        
    }
}
# ---------------------------------



# ---------------------------------
identify_gaps <- function (z, gap_tsh) {
    
    # return periods with no data (gaps) in a form of data.frame with beggining
    # and end of the period
    # Inputs:  z - zoo time series
    #          gap_tsh - interval threshold in seconds to indenrifwhich is 
    #          assumed to be a gap
    # Outputs: data.frame with begginings and ends
    
    require(zoo)
    step <- as.numeric(index(z)[-1]) - as.numeric(index(z)[-length(z)])
    
    return(data.frame(index(z[step > gap_tsh]),
                      index(z[which(step > gap_tsh) + 1])))
}
# ---------------------------------



# ---------------------------------
identify_drops <- function(x, tsh, report = T){
    
    # identify suspicious sudden dropdowns in time series
    #
    # Inputs: x - vector with values
    #         tsh - threshold to classify change as errourneous    
    #         report - should be number of of identified values reported?
    # Outputs: id.out - indexes of values identified as dropdowns
    
    x <- as.numeric(x)
    dif <- x[-1]-x[-length(x)]
    
    id_out <- which(-dif[-length(dif)] > tsh & dif[-1] > tsh) + 1
    
    if(report==T){
        print(paste("Values Identified:", length(id_out)))    
    }
    
    return(id.out)
}
# ---------------------------------



# ---------------------------------
interpolate_neighbours <- function(z, at) {
    
    # interpolate selected elements from previous and antecendent value
    #
    # Inputs: z - vector (time series)
    #         at - vector with indexes of values which should be interpolated
    # Outputs: zoo series correpsonding to z with interpolated values
    
    t <- index(z)
    d <- coredata(z)
    df <- data.frame(d[at - 1], d[at + 1])
    
    d2 <- apply(df, 1, mean, na.rm = T)
    
    d[at] <- d2
    
    return(zoo(d, t))
}
# ---------------------------------



# ---------------------------------
interpolate_drops <- function (z, tsh) {
    
    # interpolate to get rid of sudden (errourneous) dropdowns
    #
    # Inputs: z - vector (time series)
    #         tsh - threshold to classify change as errourneous
    # Outputs: id.out - indexes of values identified as dropdowns
    
    id <- identify_drops(z, tsh)
    
    if (length(id) > 0) {
        z <- interpolate_neighbours(z, id)    
    }
    
    return(z)
}
# ---------------------------------
