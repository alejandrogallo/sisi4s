#ifndef _ANGULAR_MOMENTUM_HEADER_ALE_YEAH
#define _ANGULAR_MOMENTUM_HEADER_ALE_YEAH

#include <vector>
#include <string>

namespace am {
  enum AngularMomentum { S = 1, P = 3 , D = 5 , F = 7
                       , G = 9, H = 11, I = 13, K = 15
                       };
  std::vector<AngularMomentum> all();
  size_t toInt (const AngularMomentum &);
  AngularMomentum fromString(const std::string &);
}

#endif
