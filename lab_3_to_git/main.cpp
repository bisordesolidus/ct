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
    values[0] = 1000*exp(-(dx*dy + dy*dy) / 0.02);
    values[1] = 1000*exp(-(dx*dy + dy*dy) / 0.02);
  }
};

class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0;
  }
};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    auto xx = x[0]-0.5;
    auto yy = x[1]-0.5;
    return yy+xx*(xx+pow(yy, 0.25))< DOLFIN_EPS;
  }
};

int main()
{
  auto mesh = std::make_shared<Mesh>(
    UnitSquareMesh::create({{32, 32}}, CellType::Type::triangle));
  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Constant>(0.0);
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
