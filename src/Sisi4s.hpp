#ifndef SISI4S_DEFINED
#  define SISI4S_DEFINED

#  include <util/Time.hpp>
#  include <util/Documentation.hpp>

#  include <util/CTF.hpp>
#  include <Options.hpp>

namespace sisi4s {
class Sisi4s {
public:
  void run();
  void dryRun();

  // static properties, accessible from everywhere
  static CTF::World *world;
  static Options *options;

  void printBanner();
  void printStatistics();
  void listHosts();
};
} // namespace sisi4s

#endif

/*!
 *
 * \mainpage
 * \section intro Introduction
 *
 * Coupled Cluster For Solids (SISI4S) is is a parallel quantum
 * chemistry package built on the Cyclops Tensor Framework which
 * provides high-performance structured tensor operations. SISI4S is
 * primarily focused on iterative CC methods such as CCD, CCSD, and
 * CCSD(T) for periodic systems in a plane-wave basis. The code is
 * currently interfaced to the Vienna Ab-inition Simulation Package.
 *
 * The software is available on https://github.com/alejandrogallo/sisi4s
 * and maybe obtained via the command
 *
 *
 * ~~~
 *   git clone git@github.com:alejandrogallo/sisi4s.git
 * ~~~
 *
 *
 */
