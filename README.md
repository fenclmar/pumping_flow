# pumping_flow
Some R functions to infere flow data from sewer pumping stations. Primarily designed for sewer inflow/infiltration analysis. The functions were tested on virtual data and four pumping stations in Taarnby (DK) catchment.

There are four set of functions, fun_pumping_flows.R, fun_preprocess_pumps.R, fun_Taarnby.R, and fun_virtual.R.
For the time being functions are not designed to be a standard pacakage.

## fun_pumping_flows.R
Functions for estimating inflow form current and pump sump data. Different estimation strategies are implemented as well as several optimization procedures to identify pump rules, pump capacity, etc.

## fun_preprocess_pumps.R
Functions to treat irregular times steps, missing values and outliers in real pump sump data. Functions are designed specifically for Taarnby data set.

## fun_Taarnby.R
Function for converting pump sump levels into volumes. Desingend specificaly for Taarnby pump sumps.

## fun_virtual_II.R
Functions to support realistic generating of virtual inflow/infiltration patterns. Functions are designed to work with sewage pattern generator (SPG) https://github.com/scheidan/SPG

## muskingum_simulation.cpp
Muskingum routing in C used for transformation of rainfall to runoff when generating fast and slow I/I. Coded and kindly provided by Morten Grum, WaterZerv.
