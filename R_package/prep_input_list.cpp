#include <Rcpp.h>
using namespace Rcpp;

//  Function assumes that the func_tab column names are equal (and identically ordered) to the abun_tab rownames

// [[Rcpp::export]]
List prep_input_list(NumericMatrix abun_tab, NumericMatrix func_tab) {
  
  int num_samples = abun_tab.ncol();
  int num_func = func_tab.nrow();
  
  int num_instances = num_func * num_samples;
  
  List all_sample_taxa_abun(num_instances);
  
  int instance_num = 0;
  
  for (int i = 0; i < num_func; ++i) {
    
    NumericVector func_row = func_tab(i, _); 
    LogicalVector func_row_nonzero = func_row > 0;
    
    for (int j = 0; j < num_samples; ++j) {
      
      NumericVector sample_taxa_abun = abun_tab(_, j);
      sample_taxa_abun = sample_taxa_abun[func_row_nonzero];
      LogicalVector sample_taxa_abun_nonzero_flag = sample_taxa_abun > 0;
      all_sample_taxa_abun[instance_num] = sample_taxa_abun[sample_taxa_abun_nonzero_flag];
      
      instance_num += 1;
      
    }
    
  }
  
  return all_sample_taxa_abun;
  
}
