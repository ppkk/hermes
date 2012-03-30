#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat, Hermes::Hermes1DFunction<double>* epsilon);
};
