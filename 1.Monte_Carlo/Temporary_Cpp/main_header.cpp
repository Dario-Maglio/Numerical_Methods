#include <iostream>
#include "class_fourvector.h"

int main() {
  FourVector v1(32., 40., 30., 500.);
  FourVector v2(2., 47., 60., 600.);
  Particle electron(-1, v1);
  Particle positron(+1, v2);
  std::cout<<"Total momentum: "<<(v1+v2).pt()<<std::endl;
  std::cout<<"Mass: "<<(v1+v2).mass()<<std::endl;
  electron.print();
  positron.print();
};
