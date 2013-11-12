
#include <boost/graph/detail/vebl_tree_iterators.hpp>


#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <cstdio>

using namespace boost::graph::detail;


template <std::size_t Arity>
void print_vebl_depth_records(const vebl_depth_records& rec) {
  std::vector< std::string > out_buf;
  std::cout << "    T    B    D\n";
  for(std::size_t i = 0; i < rec.T.size(); ++i) {
    out_buf.push_back(std::string(rec.T.size(),' '));
    std::cout << std::setw(5) << rec.T[i] 
              << std::setw(5) << rec.B[i] 
              << std::setw(5) << rec.D[i] << std::endl;
  };
  std::cout << std::endl;
  
  for(std::size_t i = 0; i < rec.T.size(); ++i) {
    std::size_t top_d    = get_bfl_index_depth<Arity>(rec.T[i]);
    std::size_t bottom_d = get_bfl_index_depth<Arity>(rec.B[i]);
    
    out_buf[bottom_d-1][i] = '+';
    for(std::size_t j = i+1; j < i + bottom_d; ++j)
      out_buf[bottom_d-1][j] = '-';
    
    if(i > 0)
      out_buf[top_d-1][i-top_d] = '+';
    for(std::size_t j = i-top_d+1; j < i; ++j)
      out_buf[top_d-1][j] = '-';
    
  };
  
  std::cout << "vEBL diagram:\n";
  for(std::size_t i = 0; i < out_buf.size(); ++i) {
    std::cout << out_buf[i] << std::endl;
  };
  std::cout << std::endl;
  
};

template <typename ValueType>
void print_vebl_values(const std::vector< ValueType >& values) {
  std::cout << "vEBL values:\n";
  for(std::size_t j = 0; j < values.size(); ++j)
    std::cout << std::setw(5) << values[j];
  std::cout << "\n" << std::endl;
};


int main(int argc, char** argv) {
  
  int count = 5;
  if(argc > 1) {
    count = std::atoi(argv[1]);
  };
  
  vebl_depth_records drec;
  drec.T.push_back(0); drec.B.push_back(1); drec.D.push_back(0);
  std::vector< std::size_t > values;
  values.push_back(0);
  
  print_vebl_depth_records<2>(drec);
  print_vebl_values(values);
  
  for(int i = 0; i < count; ++i) {
//     extend_vebl_depth_records<2>(drec);
    extend_vebl_storage<2>(values, drec);
    
    print_vebl_depth_records<2>(drec);
    
    std::size_t subtotal_nodes = s_treesize<2>(drec.T.size()-2);
    std::size_t total_nodes = s_treesize<2>(drec.T.size()-1);
    
    for(std::size_t j = subtotal_nodes; j < total_nodes; ++j)
      values[convert_bfl_to_vebl<2>(j,drec)] = j;
    print_vebl_values(values);
    
    std::cout << "vEBL indices:\n";
    for(std::size_t j = 0; j < total_nodes; ++j)
      std::cout << std::setw(5) << convert_bfl_to_vebl<2>(j,drec);
    std::cout << "\n" << std::endl;
  };
  
  
  
  return 0;
};








