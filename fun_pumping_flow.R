# ===========================
# Pump sump inflow estimation:
# ===========================

# Functions for estimating flow from pumping data having pump sump volumes and
# current. Exact time of pump starts and ends are not known due to raw sampling
# frequency.

# by Martin Fencl
# martinfencl@seznam.cz


get_pump_state <- function(V, current) {
    
    # function to classify pump states and returning flow table in a form
    # suitable for further processing
    # Arguments:
    # V - xts time series of PS volumes
    # current - xts time series of PS currents
    # Output:
    # table with pump states (vector with 'A', 'B', 'C', or 'D'), volumes at the
    # begginig (V0) and end of time step (v1) and duration of pump state (dt)
    # Details:
    # Pump state A is assigned to time intervals where pumping is off for the 
    # whole duration of the interval, B to the interval during which pumping starts
    # (i.e. pump is partly off and partly on), C to the intervals with pumping
    # within the whole interval, and D to the intervals during which pumping ends
    # (i.e. pump is partly on and partly off).

    # pumping on/off at time step i and (i - 1)
    on <- current > 0
    on_prev <- lag(on, k = 1)
    names(on) <- 'on'
    names(on_prev) <- 'on_prev'
    
    # Get time intervals and V(i) and V(i-1) volumes
    dt <- as.numeric(index(on)[-1]) - as.numeric(index(on)[-length(on)]) # intervals between observations
    dt <- c(NA, dt)
    v1 <- as.numeric(V)
    v0 <- as.numeric(lag(V, k = 1))
    Qin <- rep(NA, length(v1))
    tp <- dt * as.numeric(on)
    
    # classify pump states 
    pump_state <- rep(NA, length(on))
    pump_state[ which((on == F) & (on_prev == F)) ] <- 'A' # state A
    pump_state[ which((on == T) & (on_prev == F)) ] <- 'B' # state B
    pump_state[ which((on == T) & (on_prev == T)) ] <- 'C' # state C
    pump_state[ which((on == F) & (on_prev == T)) ] <- 'D' # state D
    
    ###    
    p_tab <- data.frame(pump_state, v0, v1, dt)
    colnames(p_tab) <- c('pump_state', 'V0', 'V1', 'dt')
    p_tab$pump_state <- as.character(p_tab$pump_state)
    
    return(p_tab)
    
}  


estim_next_PS_state <- function (v0, state, Qa, Qb, dt, Qp, vmin, vmax) {
    
    # function to estimate next volume in a pump sump and pump state (on/off)
    # Arguments:
    # v0 - volume (i-1) in the pump sump [l]
    # state - state is 'on' or 'off'
    # dt - time step
    # Qa - average inflow during no-pumping state
    # Qb - average iflow during pumping state
    # Qp - pump rate [l/s]
    # vmin - volume threshold for switching pum off [l]
    # vmax - volume threshold for switching pum on [l]
    
    
    # returns estimated volume and pump state at time (i)
    # (on/off)
    
    if (is.na(v0)) {
        return (c('v1' = NA, 'state'= NA, 'ton' = NA, 'toff' = NA,
                  'off_on' = NA, 'on_off' = NA))
    }
    
    
    vi <- v0    # initial flow
    t_on <- 0   # duration of pumping (pump being on)
    t_off <- 0  # duration of no-pumping (pump being off)
    t_rest <- dt - (t_off + t_on)
    if (vi > vmax) {state <- 'on'}
    if (vi < vmin) {state <- 'off'}
    t_fill <- NA
    t_empt <- NA
    on_off <- 0
    off_on <- 0
    
    while (t_rest > 0) {
        
        if (state == 'off') {
            
            # estimate filling time
            t_fill <- (vmax - vi) / Qa
            
            if (t_rest <= t_fill) {
                t_off <- t_off + t_rest
                vi <- vi + t_left * Qin
                t_rest <- 0
            } else {
                t_off <- t_off + t_fill
                t_rest <- t_rest - t_fill
                state <- 'on'
                vi <- vmax
                on_off <- on_off + 1
            }
        }    
        
        if (state == 'on') {
            if (Qin < Qp) {
                # estimate emptying time (set negatives to zero)
                t_empt <- (vi - vmin) / (Qp - Qb)
                
                if (t_rest <= t_empt) {
                    t_on <- t_on + t_rest
                    t_rest <- 0
                } else {
                    t_on <- t_on + t_empt
                    t_rest <- t_rest - t_empt
                    state <- 'off'
                    vi <- vmin
                    off_on <- off_on + 1
                }
            } else {
                t_on <- t_left
                #vi <- vi + t_left * Qin
                t_left <- 0
            }
        }
        
        # print(c('filling time' = t_fill, 'emptying time' = t_empt))
    }            
    
    v1 <- v0 + Qa * t_off + (Qb - Qp) * t_on
    
    
    return (list('v1' = v1, 'state'= state, 'ton' = t_on, 'toff' = t_off,
                 'off_on' = off_on, 'on_off' = on_off))
    
}



# ------------------
# Pump sump flows I
# ------------------

# The inflow is estiamted directly from volume and current data assuming constant
# inflow within the time step:

get_inflow_A <- function(v0, v1, dt, ...) {
    
    # get inflow to the pump sump during state A (no pumping -> no pumping)
    # Arguments:
    # v0 - volume (i-1) in the pump sump [l]
    # v1 - volume (i) in the pump sump [l]
    # dt - time difference between v1 and v0 [s]
    # Output:
    # vector (time series) with inflows [l/s]
    
    return( (v1 - v0)/dt )
}


get_inflow_B <- function(v0, v1, dt, Qp, vmax, ...) {
   
    # get inflow to the pump sump during state B (no pumping -> pumping)
    # Arguments:
    # v0 - volume (i-1) in the pump sump [l]
    # v1 - volume (i) in the pump sump [l]
    # dt - time difference between v1 and v0 [s]
    # Qp - pump rate [l/s]
    # vmax - volume threshold for switching pum on [l]
    # Output:
    # data.frame with inflows and pumping times for both solutions of
    # quadr. equation
    
    a <- dt
    b <- v0 - v1 - Qp * dt
    c <- Qp * (vmax - v0)
    
    # two solutions
    Q1 <- (-b + sqrt (b^2 - 4 * a * c)) / (2 * a)
    Q2 <- (-b - sqrt (b^2 - 4 * a * c)) / (2 * a)
    tp1 <- round((vmax - v1) / (Qp - Q1), 4)
    tp2 <- round((vmax - v1) / (Qp - Q2), 4)
    
    Q_tab <- data.frame (Q1, Q2, tp1, tp2)
    
    return( Q_tab )
}


get_inflow_B2 <- function(v0, v1, Q0, dt, Qp, vmax, ...) {

    
}



get_inflow_C <- function(v0, v1, dt, Qp, ...) {
    
    # get inflow to the pump sump during state C (pumping -> pumping)
    # Arguments:
    # v0 - volume (i-1) in the pump sump [l]
    # v1 - volume (i) in the pump sump [l]
    # dt - time difference between v1 and v0 [s]
    # Qp - pump rate [l/s]
    # Output:
    # vector (time series) with inflows [l/s]

    return( (v1 - v0) / dt + Qp )
}


get_inflow_D <- function(v0, v1, dt, Qp, vmin, ...) {
    
    # get inflow to the pump sump during state D (no pumping -> pumping)
    # Arguments:
    # v0 - volume (i-1) in the pump sump [l]
    # v1 - volume (i) in the pump sump [l]
    # dt - time difference between v1 and v0 [s]
    # Qp - pump rate [l/s]
    # vmin - volume threshold for switching pum off [l]
    # Output:
    # data.frame with inflows and pumping times for both solutions of
    # quadr. equation
    
    a <- dt
    b <- v0 - v1 - Qp * dt
    c <- Qp * (v1 - vmin)

    # two solutions
    Q1 <- (-b + sqrt (b^2 - 4 * a * c)) / (2 * a)
    Q2 <- (-b - sqrt (b^2 - 4 * a * c)) / (2 * a)
    
    
    tp1 <- round((v0 - vmin) / (Qp - Q1), 4)
    tp2 <- round((v0 - vmin) / (Qp - Q2), 4)
    
    Q_tab <- data.frame (Q1, Q2, tp1, tp2)
    
    return( Q_tab )
}



# Functions for identification of most likely solution of get_inflow_B resp. D


select_closest_solution <- function(mtx, ref) {
    
    # function to select from matrix of solutions for each row the one closest
    # to the expected value
    #
    # mtx: matrix with solutions (each column is one set of solutions)
    # ref: reference values to which mtx is compared 
    # returns vecotr with closest solutions
    
    id <- identify_closest_match(mtx, x = ref)
    
    for (i in 1:ncol(mtx)) {
        mtx[ ,i] <- mtx[ ,i] * (id == i)
    }
    
    y <- apply(mtx, 1, sum, na.rm = T)
    y[is.na(id)] <- NA
    return(y)
} 


identify_closest_match <- function(mtx, x) {
    
    # compare column of matrix with vector x (nrow(mtx) == length(x)) and
    # select the value closest to x. return vector with column numbers.
    
    if (length(x) != nrow(mtx)) {stop ('x length does not match number of mtx rows!')}
    
    difmtx <- abs(mtx - x)
    id <- numeric(length = length(x))
    id[] <- NA
    
    for (i in 1 : nrow(mtx)) {
        id[i] <- which(difmtx[i, ] == min(difmtx[i, ], na.rm = T))[1]
    }
    
    return(id)
}


get_PS_inflow <- function(xts_V, xts_cur, Qp, vmin, vmax, qref) {
    # function calculate inflow into pump sump based on volume and current data
    # Arguments:
    # xts_V - xts time series with pump sump volumes
    # xts_cur - xts time series with current data
    # dt - time difference between v1 and v0 [s]
    # Qp - pump rate [l/s]
    # vmin - volume threshold for switching pum off [l]
    # vmax - volume threshold for switching pum on [l]
    # qref - vector with 'reference' flows used during states B and D for
    #        selecting one of the solutions of quadr. eqation

    require(xts)
    pump_state <- get_pump_state (xts_V, xts_cur)
    
    idA <- which(pump_state$pump_state == 'A')
    idB <- which(pump_state$pump_state == 'B')
    idC <- which(pump_state$pump_state == 'C')
    idD <- which(pump_state$pump_state == 'D')
    
    qinA <- get_inflow_A (pump_state$V0[idA], pump_state$V1[idA],
                          pump_state$dt[idA]) 
    qinB <- get_inflow_B (pump_state$V0[idB], pump_state$V1[idB],
                          pump_state$dt[idB], Qp, vmax)
    qinC <- get_inflow_C (pump_state$V0[idC], pump_state$V1[idC],
                          pump_state$dt[idC], Qp)
    qinD <- get_inflow_D (pump_state$V0[idD], pump_state$V1[idD],
                          pump_state$dt[idD], Qp, vmin)
    
    qinB <-  select_closest_solution(qinB[, 1:2], coredata(qref)[idB])
    qinD <-  select_closest_solution(qinD[, 1:2], coredata(qref)[idD])
    
    q_est <- xts_V
    q_est[] <- NA
    q_est[idA, ] <- qinA 
    q_est[idB, ] <- qinB 
    q_est[idC, ] <- qinC 
    q_est[idD, ] <- qinD 
    q_est <- lag(q_est, k = -1)
    
    return(q_est)

}




# ------------------
# Pump sump flows II
# ------------------

# The pumping times needed for inflow estiation are estiamted from expected
# inflow, PS volumes, PS curent and pump rules. Estimation is proceed both 
# forward (from data_(i-1)) and backward (from data_(i))


tp_from_inflow <- function (Qin, xts_V, xts_cur, Qp, vmax, vmin) {
    
    # Function to calculate time of pumping based on expected inflow, PS volumes,
    # PS curent and pump rules
    # Arguments:
    # Qin - vector or time series of expected inflows
    # xts_V - xts time series with PS volumes
    # xts_V - xts time series with PS current
    # Qp - pump rate
    # vmax - volume threshold for switching the pump on
    # vmin - volume threshold for switching the pump off
    #
    # Output: data frame with two estimates of pumping times (forward and 
    #         backward) and other usefull info such as pump states, volumes, etc.
    
    if( !(length(Qin) == length(xts_V) & length(Qin) == length(xts_cur)) ) {
        stop ('Qin, xts_V and xts_cur do not have the same length!')
    }
    
    
    pump_state <- get_pump_state (xts_V, xts_cur)
    idA <- which(pump_state$pump_state == 'A')
    idB <- which(pump_state$pump_state == 'B')
    idC <- which(pump_state$pump_state == 'C')
    idD <- which(pump_state$pump_state == 'D')
    
    dt <- pump_state$dt
    
    tp_tab <- pump_state
    tp_tab$tp1 <- NA
    tp_tab$tp2 <- NA
    
    
    tp_tab$tp1[idA] <- 0
    tp_tab$tp2[idA] <- 0
    tp_tab$tp1[idB] <- dt[idB] - (vmax - tp_tab$V0[idB]) / Qin[idB]
    tp_tab$tp2[idB] <- (vmax - tp_tab$V1[idB]) / (Qp - Qin[idB])
    tp_tab$tp1[idC] <- dt[idC]
    tp_tab$tp2[idC] <- dt[idC]
    tp_tab$tp1[idD] <- (tp_tab$V0[idD] - vmin) / (Qp - Qin[idD])
    tp_tab$tp2[idD] <- dt[idD] - (tp_tab$V1[idD] - vmin) / Qin[idD]
    
    return(tp_tab)
}





# ===========================
# PS model 2:
# ===========================
# optim functions




# ===========================
# PS model fitting:
# ===========================
# optim functions

eval_2D_param_space <- function (f, p1_lim, p2_lim, step = NULL, ...) {
    
    # function to evaluate 2D parameter space (as for optimize function)
    # Arguments:
    # f - name of function (p1, p2, ...) to be evaluated
    # p1_lim - lower and upper bound of p1
    # p1_lim - lower and upper bound of p2
    # step - integer or two elemt vector step at which parameter values between
    #        lower and upper limits are evaluated
    # ... - other arguments required by function f
    
    if (is.null(step)) {
        step1 <- (p1_lim[2] - p1_lim[1]) / 10
        step2 <- (p2_lim[2] - p2_lim[1]) / 10
    } else if (length(step) == 1) {
        step1 <- step
        step2 <- step
    } else {
        step1 <- step[1]
        step2 <- step[2]
        if (length(step) > 2) {
            warning('Only two first elements of step vector are used!')
        } 
    }
    
    p1 <- seq(p1_lim[1], p1_lim[2], step1)
    p2 <- seq(p2_lim[1], p2_lim[2], step2)
    print(p1)
    print(p2)
    
    cost <- matrix(NA, length(p1), length(p2))
    
    for (i in 1 : nrow(cost)) {
        for (j in 1 : ncol(cost)) {
            cost[i, j] <- f (p = c(p1[i], p2[j]), ...)        
        }
    }    
    
    return (list('p1'= p1, 'p2' = p2, 'cost' = cost))
}


# ------------

fit_pumps <- function(ts_pump, p1.lim, p2.lim, p3.lim)
    # function to fit power law k-R model: R = alpha*(k - kw)^beta;
    # where kw = specific wet antenna attenuation
    
    # ts_pump - time series with current (first column) and volume (second column)
    # ref - vector with reference rain rates
    # p1.lim - Qp - c(min, max, ini))
    # p2.lim - Vmin - c(min, max, ini))
    # p3.lim - Vmax - c(min, max, ini))
    # p4.lim - T,F discriminat (should discriminant be added (T)or subtracted (F)?
    
{    
    
    #p.ini <- c(mean(p1.lim), mean(p2.lim)) #
    p.ini <- c(p1.lim[3], p2.lim[3], p3.lim[3])
    
    p <- optim(par = p.ini, fn = min_fun,
               lower= c(p1.lim[1], p2.lim[1], p3.lim[1]),
               upper= c(p1.lim[2], p2.lim[2], p3.lim[2]),
               method= "L-BFGS-B",
               obs = ts_pump)$par
    
    # p <- optim(par = p.ini, fn = min_fun,
    #            method= "BFGS",
    #            obs = ts_pump)$par
    
    return(p)
}

# ------------

fit_pumps_DE <- function(ts_pump, p1.lim, p2.lim, p3.lim)
    #function to fit power law k-R model: R = alpha*(k - kw)^beta;
    # where kw = specific wet antenna attenuation
    
    # ts_pump - time series with current (first column) and volume (second column)
    # ref - vector with reference rain rates
    # p1.lim - Qp - c(min, max, ini))
    # p2.lim - Vmin - c(min, max, ini))
    # p3.lim - Vmax - c(min, max, ini))
    # p4.lim - T,F discriminat (should discriminant be added (T)or subtracted (F)?
    
{    
    
    
    require(DEoptim)

    p <- DEoptim(fn = min_fun,
                     obs = ts_pump,
                     lower= c(p2.lim[1], p3.lim[1]),        
                     upper= c(p2.lim[2], p3.lim[2]),
                 control = DEoptim.control(itermax = 50))
    

    
    return(p)
}

# ------------

min_fun <- function(pars, obs)
    #function, which is minimized in fitting

{
    res <- get_inflow_volumes_v2(obs[ ,1], obs[ ,2], 9, pars[2], pars[3])

    cost_val <- cost_fun(res$Qin)
    return(cost_val)

}

# ------------

cost_fun <- function(x){
    # const function to minimize
    nas.ratio <- length(which(is.na(x))) / length(x)
    mean(abs(x - 2), na.rm = T) + 0 * nas.ratio
}


# ------------------------------
# Time series analysis
# ------------------------------


# ------------
deriv_xts <- function(x) {
    
    dval <- x - lag(x, k = 1)
    dt <- as.numeric(c(NA, index(x[-1]) - index(x[-length(x)]) ))
    
    
    return(dval/dt)
}

# ------------

get_dt <- function(x) {
    
    tim <- time(x)
    n <- length(tim)
    
    dif <- as.numeric(difftime(tim[-1], tim[-n], units = 'secs'))
    dif <- c(NA, dif)
    
    return(dif)
    
}

# ------------

d_vec <- function(x){
    dx <- x - (c(NA, x)[1 : length(x)])
    return(dx)
}

# ------------

as_naxts <- function(x) {
    
    y <- x[, 1]
    y[, 1] <- NA
    
    return(y)
    
}

# ------------
# Experimental (mot so well tested):

get_volume_peaks <- function(x){
    
    dx <- x - lag(x, k = 1)
    dx2 <- lag(dx, k = 1)
    
    idmin <- which((dx >= 0) & (dx2 < 0))
    idmax <- which((dx <= 0) & (dx2 > 0))
    
    peaks <- x
    peaks[ ,1] <- NA
    peaks[dx > 0] <- 1 # filling
    peaks[dx < 0] <- -1 # emptying
    peaks[idmin] <- -2 # minimal volume (before filling starts)
    peaks[idmax] <- 2 # maximal volume (before emptying starts)
    peaks <- lag(peaks, k = -1)
    
    return(merge(x, peaks))
    
}






