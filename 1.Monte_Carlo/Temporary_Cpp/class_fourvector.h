#include <math.h>

class FourVector {
public:
  FourVector() {}
  FourVector(float px, float py, float pz, float E) :
    px_(px), py_(py), pz_(pz), E_(E) {}
  float pt() {return sqrt(px_*px_ + py_*py_ + pz_*pz_);}
  float mass() {return sqrt(E_*E_ - (px_*px_ + py_*py_ + pz_*pz_));}
  FourVector operator+(const FourVector & other ) {
    return FourVector(px_ + other.px_, py_ + other.py_, pz_ + other.pz_, E_ + other.E_);
  };
private:
  float px_, py_, pz_;
  float E_;
};

class Particle : public FourVector{
public:
  Particle(float charge, const FourVector & p4):
    FourVector(p4), charge_(charge){}
  void print(){std::cout<<"Particle with mass "<<mass()
              <<" and charge "<<charge_<<std::endl;}
private:
  float charge_;
};
