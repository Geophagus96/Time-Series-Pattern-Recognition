#include <RcppArmadillo.h>
#include <math.h>
#include "dist_upattern.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
mat pattern_score(const vec &origin_data, const int &min_window, const int &max_window);

mat pattern_score(const vec &origin_data, const int &min_window, const int &max_window){
  int n = origin_data.n_elem;
  mat score_index((n-min_window),2,fill::ones);
  for (int i=min_window; i<n; i++){
    if (i<max_window){
      vec score((i-(min_window-1)));
      for(int j=min_window; j<=i; j++){
        vec slice = origin_data(span((i-j),(i-1)));
        double dist = dist_upattern(slice);
        score((j-min_window))=dist;
      }
      uvec index = find(score==score.min());
      score_index((i-min_window),1) = score.min();
      score_index((i-min_window),0) = index(0)+min_window;
    }
    else{
      vec score((max_window-min_window+1));
      for(int j=min_window; j<=max_window; j++){
        vec slice = origin_data(span((i-j),(i-1)));
        double dist = dist_upattern(slice);
        score((j-min_window))=dist;
      }
      uvec index = find(score==score.min());
      score_index((i-min_window),1) = score.min();
      score_index((i-min_window),0) = index(0)+min_window;
    }
  }
  return score_index;
}
