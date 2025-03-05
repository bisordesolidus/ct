#include <dolfin.h>
#include "Poisson.h"
#include <vector>

using namespace dolfin;


class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0];
    double dy = x[1];
    values[0] = 500 * exp(-(dx*dx + dy*dy) * 1000);
    values[1] = 2000 * exp(-(dx*dx + dy*dy) * 1000);
  }
};

class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 20 * exp(-1*(x[0]*x[0] + x[1]*x[1]))*sin(x[0]*x[0] + x[1]*x[1])*sin(x[0]*x[0] + x[1]*x[1]);
  }
};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    auto e = 2;
    auto xx = x[0]-5;
    auto yy = x[1];
    return pow(xx, 2/3) + pow(yy, 2/3) < pow(1, 2/3);
    //return pow(xx,3)*yy - pow(yy,3)*xx<0;
    //return xx*xx*sin(pow(yy, 2))-yy*yy*cos(pow(xx, 2))>0;
    //return xx < 100 and yy <100;
  }
};

int main()
{
  auto mesh = std::make_shared<Mesh>("figure.xml");
	auto V = std::make_shared<Poisson::FunctionSpace>(mesh);
  auto u0 = std::make_shared<Constant>(2.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  Function u(V);
  solve(a == L, u, bc);

  File file("poisson.pvd");
  file << u;

  return 0;
}