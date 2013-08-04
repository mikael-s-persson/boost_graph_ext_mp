
#include <vector>
#include <list>
#include <iostream>


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



int main() {
  
  std::list<int> l;
  l.push_back(42);
  l.push_back(69);
  
  std::vector<int> v;
  v.push_back(42);
  
  foo(l);
  foo(v);
  foo(&l);
  foo(&v);
  
  return 0;
};








