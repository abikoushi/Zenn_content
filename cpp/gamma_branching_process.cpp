#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame simulate_branching_cpp(double a, double b, double Tmax) {
  
  std::vector<int> id;
  std::vector<int> parent;
  std::vector<double> birth_time;
  std::vector<int> generation;
  
  id.push_back(0);
  parent.push_back(NA_INTEGER);
  birth_time.push_back(0.0);
  generation.push_back(0);
  
  int next_id = 1;
  size_t i = 0;
  
  while (i < id.size()) {
    
    double birth = birth_time[i];
    double t = birth;
    
    while (true) {
      
      t += R::rgamma(a, 1.0 / b);
      
      if (t > Tmax) {
        break;
      }
      
      id.push_back(next_id);
      parent.push_back(id[i]);
      birth_time.push_back(t);
      generation.push_back(generation[i] + 1);
      
      ++next_id;
    }
    
    ++i;
  }
  

  std::vector<size_t> order(birth_time.size());
  for (size_t j = 0; j < order.size(); ++j) {
    order[j] = j;
  }
  std::sort(
    order.begin(),
    order.end(),
    [&](size_t x, size_t y) {
      return birth_time[x] < birth_time[y];
    }
  );
  

  IntegerVector id_out(order.size());
  IntegerVector parent_out(order.size());
  NumericVector birth_out(order.size());
  IntegerVector generation_out(order.size());
  for (size_t j = 0; j < order.size(); ++j) {
    
    size_t k = order[j];
    
    id_out[j] = id[k];
    parent_out[j] = parent[k];
    birth_out[j] = birth_time[k];
    generation_out[j] = generation[k];
  }
  
  return DataFrame::create(
    _["id"] = id_out,
    _["parent"] = parent_out,
    _["birth_time"] = birth_out,
    _["generation"] = generation_out
  );
}