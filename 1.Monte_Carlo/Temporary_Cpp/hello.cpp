//Preprocessor directives
#include <iostream>
#include <map>
#define MIN(a,b) ((a<b)?(a):(b))
// The “preprocessor” is run on the code before the actual compiler.
// It replaces the effective source code seen by the compiler

//using namespace std; to avoid std each time

int fun(int x){
  // A function's signature includes the function's name and the number,
  // order and type of its formal parameters. Two overloaded functions must
  // not have the same signature.
  return MIN(x, 42);
}

inline int fun(int x, int y) {
  // Inline make sense for small functions avoiding jumps
  return MIN(x, y);
}

int main() {
  std::cout << "Hello world!" << std::endl;
  std::cout << "The Min is: " << fun(12) << " or " << fun(12, 4) << std::endl;

  float a[] = {10,2,4,5};
  typedef float * floatiterator ;
  floatiterator begin= a;
  floatiterator end= a+4;
  for(floatiterator it=begin; it!=end;it++){
    std::cout << *it << " at position " << (it-begin) << std::endl;
  }

  std::map<int,float> dict;
  dict[10]=-1.;
  dict[7]=1.13;
  std::cout << "item 2 : " << dict[2] << std::endl;
  for(std::map<int,float>::const_iterator it=dict.begin();it!=dict.end(); it++){
    std::cout << it->first << " => " << it->second << std::endl;
  }
}

// g++ -o main Hello.cpp
// ./main
