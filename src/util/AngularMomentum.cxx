#include <util/AngularMomentum.hpp>

std::vector<am::AngularMomentum> am::all() { return {S, P, D, F, G, H, I, K}; }
size_t am::toInt(const am::AngularMomentum &a) { return a; }

am::AngularMomentum am::fromString(const std::string &a) {
  if (a == "S") return S;
  if (a == "P") return P;
  if (a == "D") return D;
  if (a == "F") return F;
  if (a == "G") return G;
  if (a == "H") return H;
  if (a == "I") return I;
  if (a == "K") return K;
  throw "I don't anderstand symbol: " + a;
}
