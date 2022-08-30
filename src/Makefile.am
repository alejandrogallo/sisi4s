AUTOMAKE_OPTIONS = subdir-objects
include $(top_srcdir)/Cc3s.mk

AM_CXXFLAGS = $(CTF_CPPFLAGS) \
              -fmax-errors=1 \
              -D_POSIX_C_SOURCE=200112L     \
              -D__STDC_LIMIT_MACROS         \
              -DFTN_UNDERSCORE=1            \
              -DCC4S_VERSION=\"cc3s\" \
              "-DCC4S_DATE=\"${DATE}\""     \
              "-DCOMPILER_VERSION=\"${COMPILER_VERSION}\"" \
              "-DCC4S_MPI_VERSION=\"${CC4S_MPI_VERSION}\""


bin_PROGRAMS = cc3s

cc3s_CPPFLAGS = $(CTF_CPPFLAGS) $(YAML_CPPFLAGS) $(LIBINT_CPPFLAGS) $(EIGEN_CPPFLAGS)

cc3s_SOURCES =                                                   \
./Data.cxx                                                       \
./Cc4s.cxx                                                       \
./turbomole/MosParser.cxx                                        \
./Options.cxx                                                    \
./util/BlacsWorld.cxx                                            \
./util/Emitter.cxx                                               \
./util/Parsing.cxx                                               \
./util/ScaLapackMatrix.cxx                                       \
./util/Log.cxx                                                   \
./util/TensorIo.cxx                                              \
./util/AngularMomentum.cxx                                       \
./util/FlopsCounter.cxx                                          \
./util/BasisSet.cxx                                              \
./util/Timer.cxx                                                 \
./nwchem/MovecsParser.cxx                                        \
./nwchem/BasisParser.cxx                                         \
./mixers/LinearMixer.cxx                                         \
./mixers/DiisMixer.cxx                                           \
./mixers/Mixer.cxx                                               \
./tcc/DryTensor.cxx                                              \
./algorithms/Algorithm.cxx                                       \
./Parser.cxx                                                     \
./math/CanonicalPolyadicDecomposition.cxx                        \
./math/ComplexTensor.cxx                                         \
./math/IterativePseudoInverse.cxx                                \
./math/PseudoInverseHermitianSvd.cxx                             \
./math/RegularizedAlternatingLeastSquares.cxx                    \
./math/PseudoInverseSvd.cxx                                      \
./math/MathFunctions.cxx                                         \
./algorithms/PseudoInverse.cxx                                   \
./algorithms/SimilarityTransformedHamiltonian.cxx                \
./algorithms/sthamiltonian/wabij.cxx                             \
./algorithms/sthamiltonian/wabcijk.cxx                           \
./algorithms/ClusterSinglesDoublesAlgorithm.cxx                  \
./algorithms/CcsdPreconditioner.cxx                              \
./algorithms/ParticleHoleCoulombIntegrals.cxx                    \
./algorithms/CcsdEnergyFromCoulombIntegralsReference.cxx         \
./algorithms/BasisSetExtrapolation.cxx                           \
./algorithms/TensorGetMax.cxx                                    \
./algorithms/TensorWriter.cxx                                    \
./algorithms/ParticleHoleCoulombVertexSingularVectors.cxx        \
./algorithms/NaturalTransitionOrbitalsFromRhoAI.cxx              \
./algorithms/Mp2ImaginaryFrequencyGrid.cxx                       \
./algorithms/ComplexTensorNorm.cxx                               \
./algorithms/GenerateRandomMatrix.cxx                            \
./algorithms/TensorSlicer.cxx                                    \
./algorithms/FockMatrixFromCoulombIntegrals.cxx                  \
./algorithms/ParticleHoleCoulombVertexReader.cxx                 \
./algorithms/TensorReader.cxx                                    \
./algorithms/MoReader.cxx                                        \
./algorithms/CoulombVertexReader.cxx                             \
./algorithms/ComplexTensorReader.cxx                             \
./algorithms/FcidumpReader.cxx                                   \
./algorithms/FcidumpWriter.cxx                                   \
./algorithms/SliceCoulombVertex.cxx                              \
./algorithms/ComplexTensorWriter.cxx                             \
./algorithms/TensorUnrestricter.cxx                              \
./algorithms/StantonIntermediatesUCCSD.cxx                       \
./algorithms/ApproximateCoulombVertex.cxx                        \
./algorithms/CcsdPerturbativeTriples.cxx                         \
./algorithms/TensorNorm.cxx                                      \
./algorithms/TensorSum.cxx                                       \
./algorithms/UnrestrictedEquationOfMotionSinglesFromRpa.cxx      \
./algorithms/UrpaAmplitudesFromCoulombIntegrals.cxx              \
./algorithms/OneBodyRotation.cxx                                 \
./algorithms/HartreeFockFromCoulombIntegrals.cxx                 \
./algorithms/Delete.cxx                                          \
./algorithms/Mp2NaturalOrbitals.cxx                              \
./algorithms/MeanCorrelationHoleDepth.cxx                        \
./algorithms/FiniteSizeCorrection.cxx                            \
./algorithms/TensorAntisymmetrizer.cxx                           \
./algorithms/UccsdAmplitudesFromCoulombIntegrals.cxx             \
./algorithms/UCcsdIPEquationOfMotionDavidson.cxx                 \
./algorithms/UCcsdEAEquationOfMotionDavidson.cxx                 \
./algorithms/CcsdEnergyFromCoulombIntegrals.cxx                  \
./algorithms/CoulombIntegralsFromRotatedCoulombIntegrals.cxx     \
./algorithms/TensorReduceijij.cxx                                \
./algorithms/CcsdDiagrammaticDecomposition.cxx                   \
./algorithms/GenerateRandomTensor.cxx                            \
./algorithms/GenerateRandomComplexMatrix.cxx                     \
./algorithms/Exit.cxx

## FOR LIBINT2
cc3s_SOURCES +=                                                  \
./algorithms/HartreeFockFromGaussian.cxx                         \
./algorithms/OneBodyFromGaussian.cxx                             \
./algorithms/NuclearRepulsionEnergy.cxx                          \
./algorithms/CoulombIntegralsFromGaussian.cxx



cc3s_LDADD = $(CTF_LDFLAGS) $(YAML_LDFLAGS) $(BLAS_LIBS) -lscalapack -lquadmath $(LIBINT_LIB)



## ./algorithms/UPerturbativeTriples.cxx                            \
## ./algorithms/DoublesAmplitudesFromVertex.cxx                     \
## ./algorithms/ParticleHoleCoulombVertexFromFactors.cxx            \
## ./algorithms/ComplexTensorSum.cxx                                \
## ./algorithms/RandomOrbitals.cxx                                  \

## ./algorithms/TensorContraction.cxx                               \
## ./algorithms/CoulombVertexFromFactors.cxx                        \
## ./algorithms/ApproximateParticleHoleCoulombVertex.cxx            \

## ./algorithms/PerturbativeTriples.cxx                             \
## ./algorithms/CoulombVertexDecomposition.cxx                      \
## ./algorithms/ParticleHoleCoulombVertexDecomposition.cxx          \
## ./algorithms/CoulombIntegralsFromVertex.cxx                      \
## ./algorithms/DoublesAmplitudesDecomposition.cxx                  \
## ./algorithms/CcsdPerturbativeTriplesComplex.cxx                  \
## ./algorithms/CcsdEquationOfMotionDavidson.cxx                    \
## ./algorithms/ComplexTensorContraction.cxx                        \
## ./algorithms/ParenthesisTriples.cxx                              \
## ./algorithms/FromComplexTensor.cxx                               \
## ./algorithms/SingleParticleOccupancies.cxx                       \
## ./algorithms/DrccdEnergyFromCoulombIntegrals.cxx                 \
## ./algorithms/PolarizabilityFromCoulombVertex.cxx                 \
## ./algorithms/PartitionTensor.cxx                                 \
## ./algorithms/Mp2EnergyFromCoulombIntegrals.cxx                   \
## ./algorithms/CoulombVertexSingularVectors.cxx                    \
## ./algorithms/ClusterSinglesDoublesTriplesQuadruplesAlgorithm.cxx \
## ./algorithms/ClusterSinglesDoublesTriplesAlgorithm.cxx           \
