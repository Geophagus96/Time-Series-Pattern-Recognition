#pragma once
#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double dist_upattern(const vec &input);

//[[Rcpp::export]]
double dist_upattern2(const vec &input);

//[[Rcpp::export]]
double dist_lpattern(const vec &input);

//[[Rcpp::export]]
double stand_dist_2series(const vec &input1, const vec &input2);

double dist_upattern2(const vec &input){
  double dist;
  if (input.has_inf()){
    dist = 1;
  }
  else{
    vec input_vec = cumprod((1+input));
    int n = input.n_elem;
    double input_max = input_vec.max();
    double input_min = input_vec.min();
    if (input_max == input_min){
      dist = 1;
    }
    else{
      vec transform_input = (1/(input_max-input_min))*input_vec-(input_min/(input_max-input_min));
      double a = 1/(pow(((n-1)/2),2));
      vec base = linspace(-((n-1)/2), ((n-1)/2),n);
      vec benchmark = a*pow(base,2);
      vec yk = transform_input(span(0,(n-2)));
      vec ykp = transform_input(span(1,(n-1)));
      vec xk = benchmark(span(0,(n-2)));
      vec xkp = benchmark(span(1,(n-1)));
      vec diff = (ykp-yk)-(xkp-xk);
      dist = as_scalar(sum(pow(diff,2))/(n-1));
    }
  }
  return dist;
}

double dist_upattern(const vec &input){
  double dist;
  if (input.has_inf()){
    dist = 1;
  }
  else{
    vec input_vec = cumprod((1+input));
    int n = input.n_elem;
    double input_max = input_vec.max();
    double input_min = input_vec.min();
    if (input_max == input_min){
      dist = 1;
    }
    else{
      vec transform_input = (1/(input_max-input_min))*input_vec-(input_min/(input_max-input_min));
      double a = 1/((n-1)/2);
      vec base = linspace(-((n-1)/2),((n-1)/2),n);
      vec benchmark = a*abs(base);
      vec yk = transform_input(span(0,(n-2)));
      vec ykp = transform_input(span(1,(n-1)));
      vec xk = benchmark(span(0,(n-2)));
      vec xkp = benchmark(span(1,(n-1)));
      vec diff = (ykp-yk)-(xkp-xk);
      dist = as_scalar(sum(pow(diff,2))/(n-1));
    }
  }
  return dist;
}

double dist_lpattern(const vec &input){
  double dist;
  if (input.has_inf()){
    dist = 1;
  }
  else{
    vec input_vec = cumprod((1+input));
    int n = input.n_elem;
    double input_max = input_vec.max();
    double input_min = input_vec.min();
    if (input_max == input_min){
      dist = 1;
    }
    else{
      vec transform_input = (1/(input_max-input_min))*input_vec-(input_min/(input_max-input_min));
      double a = 1/pow((n-1),2);
      vec base = linspace(0,(n-1),n);
      vec benchmark = a*pow(base,2);
      vec yk = transform_input(span(0,(n-2)));
      vec ykp = transform_input(span(1,(n-1)));
      vec xk = benchmark(span(0,(n-2)));
      vec xkp = benchmark(span(1,(n-1)));
      vec diff = (ykp-yk)-(xkp-xk);
      dist = as_scalar(sum(pow(diff,2))/(n-1));
    }
  }
  return dist;
}

double stand_dist_2series(const vec &input, const vec &compare){
  double dist;
  if ((input.has_inf())or(input.has_inf())){
    dist = 1;
  }
  else{
    vec input_vec = cumprod((1+input));
    int n = input.n_elem;
    double input_max = input_vec.max();
    double input_min = input_vec.min();
    double compare_max = compare.max();
    double compare_min = compare.min();
    if ((input_max == input_min)or(compare_max == compare_min)){
      dist = 1;
    }
    else{
      vec transform_input = (1/(input_max-input_min))*input_vec-(input_min/(input_max-input_min));
      vec yk = transform_input(span(0,(n-2)));
      vec ykp = transform_input(span(1,(n-1)));
      vec xk = compare(span(0,(n-2)));
      vec xkp = compare(span(1,(n-1)));
      vec diff = (ykp-yk)-(xkp-xk);
      dist = as_scalar(sum(pow(diff,2))/(n-1));
    }
  }
  return dist;
}
