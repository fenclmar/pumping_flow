#include <Rcpp.h>
using namespace Rcpp;

int last_index = 0;

IntegerVector rain_series_index;
NumericVector rain_series;

double rainValue( double t)
{
  //std::cout << "Calling rain(" << std::setprecision(32) << t << ")" << std::endl;
  if(t < rain_series_index[0] || t > rain_series_index[rain_series_index.size()-1]){
    //std::cout << "o";
    return 0.0;
  } else {
    //std::cout << "x";
    double t_upper = rain_series_index[last_index];
    double t_lower = rain_series_index[last_index-1];
    while(t > t_upper && last_index < (rain_series_index.size()-2)){
      //std::cout << "+";
      last_index = last_index + 1;
      t_upper = rain_series_index[last_index];
      t_lower = rain_series_index[last_index-1];
    }
    while(t <= t_lower && last_index > 1){
      //std::cout << "-";
      last_index = last_index - 1;
      t_upper = rain_series_index[last_index];
      t_lower = rain_series_index[last_index-1];
    }
    
    return rain_series[last_index];
  }
}

// [[Rcpp::export]]
NumericVector muskingum_reservoirs_simulation_rcpp(
    int res_count, double retention_time, IntegerVector rain_index, NumericVector rain, IntegerVector flow_index) {

  rain_series_index = rain_index;
  rain_series = rain;
  
  NumericVector flow_out(flow_index.size());
  NumericVector c(3);
  
  if(flow_index[0] < rain_index[0]){
    std::cout << "Rain should start before flow." << std::endl;
  }
  
  int run_from_time;
  int r = 0;
  int f = 0;
  if(flow_index[0] == rain_index[0]){
    run_from_time = rain_index[0];
    r = 1;
    f = 1;
  } else if(flow_index[0] > rain_index[0]){
    run_from_time = rain_index[0];
    r = 1;
    f = 0;
  } else {
    run_from_time = flow_index[0];
    r = 0;
    f = 1;
  }
  
  // Run to time for first step ...
  int run_to_time = std::min(rain_index[r], flow_index[f]);
  if(run_to_time == rain_index[r])
  {
    r++;
  }
  if(run_to_time == flow_index[f])
  {
    f++;
  }

  NumericMatrix flowValue(res_count,2);
  for(int res=0; res < res_count; res = res + 1){
    flowValue(res,0) = 0.001;
    flowValue(res,1) = 0.001;
  }
  
  
  double k = (retention_time/(res_count-1));
  double x = 0.0; // Muskingum X
  while((run_to_time < rain_index[rain_index.size()-1]) && (run_to_time < flow_index[flow_index.size()-1])){

    int delta_t = run_to_time - run_from_time;
    double denominator = 2.0*k*(1-x) + delta_t;
    c[0] = ( delta_t - 2.0*k*x )/denominator;
    c[1] = ( delta_t + 2.0*k*x )/denominator;
    c[2] = ( 2.0*k*(1-x) - delta_t )/denominator;
    
    flowValue(0, 1) = 
      c[0]*rainValue(run_to_time) + 
      c[1]*rainValue(run_from_time) + 
      c[2]*flowValue(0, 0);
    for(int res=1; res < res_count; res = res + 1){
      flowValue(res, 1) = 
        c[0]*flowValue(res-1, 1) + 
        c[1]*flowValue(res-1, 0) + 
        c[2]*flowValue(res, 0);
    }
    
    if(run_to_time == flow_index[f-1])
    {
      //std::cout << run_to_time << "\t" << flowValue(res_count-1, 1) << std::endl;
      flow_out[f] = flowValue(res_count-1, 1);
      if(flow_out[f] < 0.0){
        flow_out[f] = 0.0;
      }
    }

    run_from_time = run_to_time;
    run_to_time = std::min(rain_index[r], flow_index[f]);
    if(run_to_time == rain_index[r])
    {
      r++;
    }
    if(run_to_time == flow_index[f])
    {
      f++;
    }
    
    for(int res=0; res < res_count; res = res + 1){
      flowValue(res,0) = flowValue(res,1);
      if(flowValue(res,0) < 0.0){
        flowValue(res,0) = 0.0;
      }
    }
  }

  return flow_out;
  
}

// R Trials
/*** R
  muskingum_reservoirs_simulation <- function(res_count, retention_time, xts_rain, xts_flow){
    return(muskingum_reservoirs_simulation_rcpp(res_count, retention_time,
                                            as.numeric(index(xts_rain)), 
                                            as.numeric(coredata(xts_rain)),
                                            as.numeric(index(xts_flow))))
  }

*/
