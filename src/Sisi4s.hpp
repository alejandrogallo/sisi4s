/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SISI4S_DEFINED
#define SISI4S_DEFINED

#include <util/Time.hpp>

#include <util/CTF.hpp>
#include <Options.hpp>

namespace sisi4s {
  class Sisi4s {
  public:
    Sisi4s();
    ~Sisi4s();
    void run();
    void dryRun();

    static bool isDebugged();

    // static properties, accessible from everywhere
    static CTF::World *world;
    static Options *options;

  protected:
    void printBanner();
    void printStatistics();
    void listHosts();
  };
}

#endif

/*! \mainpage 
 * \section intro Introduction
 * 
 * Coupled Cluster For Solids (SISI4S) is is a parallel quantum
 * chemistry package built on the Cyclops Tensor Framework which
 * provides high-performance structured tensor operations. SISI4S is
 * primarily focused on iterative CC methods such as CCD, CCSD, and
 * CCSD(T) for periodic systems in a plane-wave basis. The code is
 * currently interfaced to the Vienna Ab-inition Simulation Package.
 * 
 * The software is available on
 * https://hulldamage.org/sisi4s.org/doku.php and maybe obtained
 * via the command 
 *
 * git clone git@gitlab.com:grueneis/sisi4s.git
 *
 * For information or access to the repository please contact <mailto:andreas.grueneis@tuwien.ac.at>
 *
 * SISI4S depends on the CTF and requires MPI to be built as the main parallel execution and communication mechanism.
 *
 * \section algorithms Algorithms
 *
 *The main functionality is based on different algorithms.
 */