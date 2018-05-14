/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SharedPointer.hpp>
#include <string>
#include <ctf.hpp>
#include <tcc/DryTensor.hpp>

namespace cc4s {
  /**
   * \brief Provides the functionality for a finite temperature algorithm with
   * only doubles amplitudes.
   */
  class ThermalClusterDoublesAlgorithm: public Algorithm {
  public:
    ThermalClusterDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalClusterDoublesAlgorithm();
    /**
     * \brief Calculates the energy of this ThermalClusterDoubles algorithm
     */
    virtual void run();

    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();
    /**
     * \brief Returns the abbreviation of the concrete algorithm, without
     * the leading "Thermal", e.g. "Ccd" for "ThermalCcdEnergy"
     */
    virtual std::string getAbbreviation() = 0;

    static constexpr int DEFAULT_MAX_ITERATIONS = 12;

  protected:
    /**
     * \brief doubles amplitudes on the imaginary time grid
     **/
    std::vector<PTR(CTF::Tensor<real>)> Tabijn;

    /**
     * \brief The eigenenergy difference of between the state a,b and i,j.
     * This is used to propagate all states in imaginary time.
     **/
    PTR(CTF::Tensor<>) Dabij;

    /**
     * \brief Inverse temperature \f$\beta=1/k_{\rm B}T\f$, where
     * \f$k_{\rm B}T\f$ is given in the same unit as the eigenenergies
     * \f$\varepsilon_p\f$.
     **/
    real beta;


    /**
     */
    virtual void applyHamiltonian(
      CTF::Tensor<real> &T0abij,
      CTF::Tensor<real> &T1abij,
      const real DTau,
      CTF::Tensor<real> &S1abij
    ) = 0;

    std::string getCapitalizedAbbreviation();
    std::string getAmplitudeIndices(CTF::Tensor<real> &T);
    void fetchDelta(CTF::Tensor<real> &Delta);
    void thermalContraction(CTF::Tensor<real> &T);

    class ImaginaryTimeTransform {
    protected:
      ImaginaryTimeTransform(
        real DTau_
      ): DTau(DTau_) {
      }
      real DTau;
    };

    class ImaginaryTimePropagation: public ImaginaryTimeTransform {
    public:
      ImaginaryTimePropagation(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real DeltaJ, real &T) const {
        T *= std::exp(-DeltaJ*DTau);
      }
    };

    // constant V^J=C contribution to T'^J(tau_m),
    // independent of the amplitudes T^I(tau)
    // = f^(J/I) * int_0^Tau dtau exp(-DeltaJ*(Tau-tau))
    class ConvolutionC: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      ConvolutionC(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^J
      void operator ()(const real DeltaJ, real &T) const {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1 - x*(1./2 - x*(1./6 - x/24)) );
        } else {
          T *= (1-std::exp(-x))/DeltaJ;
        }
      }
    };

    // linear T^I(tau_m-1)=T0 contribution to T'^J(tau_m)
    // = f^(J/I) * int_0^Tau dtau (Tau-tau)/Tau * exp(-DeltaJ*(Tau-tau))
    class Convolution0: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution0(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real DeltaJ, real &T) {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./2 - x*(1./3 - x*(1./8 - x/30)) );
        } else {
          T *= (1-(1+x)*std::exp(-x)) / (x*DeltaJ);
        }
      }
    };

    // linear T^I(tau_m)=T1 contribution to T'^J(tau_m)
    // = f^(J/I) * int_0^Tau dtau tau/Tau * exp(-DeltaJ*(Tau-tau))
    class Convolution1: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution1(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real DeltaJ, real &T) {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./2 - x*(1./6 - x*(1./24 - x/120)) );
        } else {
          T *= (1-(1-std::exp(-x))/x) / DeltaJ;
        }
      }
    };

    // quadratic T^I1(tau_m-1)*T^I2(tau_m-1)=T0*T0 contribution to T'^J(tau_m)
    // = f^(J/I)
    // * int_0^Tau dtau (Tau-tau)/Tau * (Tau-tau)/Tau * exp(-DeltaJ*(Tau-tau))
    class Convolution00: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution00(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real DeltaJ, real &T) {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./3 - x*(1./4 - x*(1./10 - x/36)) );
        } else {
          T *= (2-((2+2*x+x*x)*std::exp(-x))) / (x*x*DeltaJ);
        }
      }
    };

    // quadratic T^I1(tau_m-1)*T^I2(tau_m)=T0*T1 contribution to T'^J(tau_m)
    // = f^(J/I)
    // * int_0^Tau dtau (Tau-tau)/Tau * tau/Tau * exp(-DeltaJ*(Tau-tau))
    class Convolution01: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution01(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real DeltaJ, real &T) {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./6 - x*(1./12 - x*(1./40 - x/180)) );
        } else {
          T *= ((x-2)+((x+2)*std::exp(-x))) / (x*x*DeltaJ);
        }
      }
    };

    // quadratic T^I1(tau_m)*T^I2(tau_m)=T1*T1 contribution to T'^J(tau_m)
    // = f^(J/I)
    // * int_0^Tau dtau tau/Tau * tau/Tau * exp(-DeltaJ*(Tau-tau))
    class Convolution11: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution11(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real DeltaJ, real &T) {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./3 - x*(1./12 - x*(1./60 - x/360)) );
        } else {
          T *= ((2-2*x+x*x)-(2*std::exp(-x)))/ (x*x*DeltaJ);
        }
      }
    };
  };
}

#endif

