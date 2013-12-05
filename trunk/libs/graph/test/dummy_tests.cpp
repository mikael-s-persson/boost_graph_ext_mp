
#include <vector>
#include <list>
#include <iostream>

#include <limits>

#include <boost/variant.hpp>


template <typename Container>
void foo(Container& cont) {
  std::cout << "General template, with container size: " << cont.size() << std::endl;
};

template <typename ValueType>
void foo(std::vector<ValueType>& cont) {
  std::cout << "Specific template, with vector size: " << cont.size() << std::endl;
};

template <typename Container>
void foo(Container* cont) {
  std::cout << "Pointer template, with pointee container size: " << cont->size() << std::endl;
};




template <typename IntType, std::size_t Base>
struct s_power_array {
  static const std::size_t max_power = (std::numeric_limits<IntType>::digits * std::numeric_limits<IntType>::radix) / Base;
  IntType values[max_power];
  s_power_array() {
    values[0] = 1;
    for(std::size_t i = 1; i < max_power; ++i)
      values[i] = values[i-1] * Base;
  };
};

template <typename IntType, std::size_t Base>
IntType s_power(std::size_t aPower) {
  static const s_power_array<IntType, Base> powers;
  return powers.values[aPower];
};




template <typename IntType, std::size_t Base>
struct s_treesize_array {
  static const std::size_t max_depth = (std::numeric_limits<IntType>::digits * std::numeric_limits<IntType>::radix) / Base;
  IntType values[max_depth];
  s_treesize_array() {
    values[0] = 1;
    for(std::size_t i = 1; i < max_depth; ++i)
      values[i] = values[i-1] * Base + 1;
  };
};

template <typename IntType, std::size_t Base>
IntType s_treesize(std::size_t aDepth) {
  static const s_treesize_array<IntType, Base> treesizes;
  return treesizes.values[aDepth];
};



int main() {
  
  std::cout << " Sizeof int = " << sizeof(int) << " and Sizeof void* = " << sizeof(void*) << std::endl
            << " Sizeof boost::variant<int> = " << sizeof(boost::variant<int>) << std::endl
            << " Sizeof boost::variant<int, void*> = " << sizeof(boost::variant<int,void*>) << std::endl
            << " Sizeof boost::variant<int, unsigned int> = " << sizeof(boost::variant<int,unsigned int>) << std::endl;
  
  std::list<int> l;
  l.push_back(42);
  l.push_back(69);
  
  std::vector<int> v;
  v.push_back(42);
  
  foo(l);
  foo(v);
  foo(&l);
  foo(&v);
  
  std::cout << "Max power 2 of std::size_t: " << s_power_array<std::size_t, 2>::max_power << std::endl;
  std::cout << "Max power 3 of std::size_t: " << s_power_array<std::size_t, 3>::max_power << std::endl;
  std::cout << "Max power 4 of std::size_t: " << s_power_array<std::size_t, 4>::max_power << std::endl;
  std::cout << "Max power 5 of std::size_t: " << s_power_array<std::size_t, 5>::max_power << std::endl;
  std::cout << "Max power 6 of std::size_t: " << s_power_array<std::size_t, 6>::max_power << std::endl;
  std::cout << "Max power 7 of std::size_t: " << s_power_array<std::size_t, 7>::max_power << std::endl;
  std::cout << "Max power 8 of std::size_t: " << s_power_array<std::size_t, 8>::max_power << std::endl;
  
  s_power_array<std::size_t, 8> powers;
  for(auto x : powers.values)
    std::cout << x << " ";
  std::cout << std::endl;
  
  s_treesize_array<std::size_t, 8> treesizes;
  for(auto x : treesizes.values)
    std::cout << x << " ";
  std::cout << std::endl;
  
  return 0;
};








