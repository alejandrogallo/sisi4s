dnl  This is an autoconf script.
dnl  To rebuild the `configure' script from this, execute the command
dnl 	autoconf
dnl  in the directory containing this script.
dnl

AC_INIT([cc3s],[0.0.0],[https://gihtub.com/alejandrogallo/cc3s])
AC_CONFIG_SRCDIR([configure.ac])
AC_PREREQ(2.69)
AC_CONFIG_MACRO_DIRS([./etc/m4])
AC_CONFIG_AUX_DIR([./etc/autotools])
AM_INIT_AUTOMAKE([foreign -Wall])


dnl -----------------------------------------------------------------------
dnl Configure switches

dnl CTF WITH OPTION
AC_ARG_WITH(ctf,
            [AS_HELP_STRING([--with-ctf],
                           [prefix for CTF includes and libraries] )],
            [LIBCTF_PATH="`readlink -f $withval`";
             LIBCTF_CPATH="`readlink -f $withval`/include";
             LIBCTF_LD_LIBRARY_PATH="`readlink -f $withval`/lib";
             BUILD_CTF=NO],
            [BUILD_CTF=YES])

dnl LIBINT WITH OPTION
AC_ARG_WITH(libint,
            [AS_HELP_STRING([--with-libint],
                            [prefix for libint includes and libraries] )],
            [LIBINT_PATH="`readlink -f $withval`";
             LIBINT_CPATH="`readlink -f $withval`/include";
             LIBINT_LD_LIBRARY_PATH="`readlink -f $withval`/lib";
             BUILD_LIBINT=NO],
            [BUILD_LIBINT=YES])

dnl EIGEN WITH OPTION
AC_ARG_WITH(libint,
            [AS_HELP_STRING([--with-libint],
                            [prefix for libint includes and libraries] )],
            [EIGEN_PATH="`readlink -f $withval`";
             EIGEN_CPATH="`readlink -f $withval`/include";
             BUILD_EIGEN=NO],
            [BUILD_EIGEN=YES])

dnl LIBYAML WITH OPTION
AC_ARG_WITH(yamlcpp,
            [AS_HELP_STRING([--with-yamlcpp],
                            [prefix for yamlcpp includes and libraries] )],
            [YAMLCPP_PATH="`readlink -f $withval`";
             YAMLCPP_CPATH="`readlink -f $withval`/include";
             YAMLCPP_LD_LIBRARY_PATH="`readlink -f $withval`/lib";
             BUILD_YAMLCPP=NO],
            [BUILD_YAMLCPP=YES])

AC_ARG_WITH([clang-check],
            [AS_HELP_STRING([--with-clang-check],
                            [Check with clang-check])],
            [clang_check=YES;
             AC_CHECK_PROG([CLANG_CHECK], [clang-check])],
            [clang_check=NO])
AM_CONDITIONAL([WITH_CLANG_CHECK], [test x${clang_check} = xYES])

AC_ARG_ENABLE([cuda],
              [AS_HELP_STRING([--enable-cuda],
                              [Build with cuda])],
                              [WITH_CUDA=yes],
                              [WITH_CUDA=no])
AC_ARG_VAR([NVCC], [Path to the nvidia cuda compiler.])
AC_ARG_VAR([CUDA_LDFLAGS], [LDFLAGS to find libraries -lcuda, -lcudart, -lcublas.])
AC_ARG_VAR([CUDA_CXXFLAGS], [CXXFLAGS to find the CUDA headers])


dnl -----------------------------------------------------------------------
AM_CONDITIONAL([WITH_BUILD_CTF], [test x${BUILD_CTF} = xYES])
AM_CONDITIONAL([WITH_BUILD_LIBINT], [test x${BUILD_LIBINT} = xYES])
AM_CONDITIONAL([WITH_BUILD_YAML], [test xYES = xYES])
AM_CONDITIONAL([WITH_BUILD_EIGEN], [test xYES = xYES])
AC_SUBST(LIBCTF_PATH)
AC_SUBST(LIBCTF_CPATH)
AC_SUBST(LIBCTF_LD_LIBRARY_PATH)

dnl -----------------------------------------------------------------------
AC_MSG_NOTICE([

              CHECKING YOUR OVERALL SYSTEM ...

               ])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

dnl Guess compilers
AC_LANG([C++])
AC_MSG_NOTICE([
               CHECKING COMPILERS...
               ])
AC_PROG_CXX
AX_MPI([], [AC_MSG_FAILURE([No MPI C++ compiler found.])])
AX_COMPILER_VENDOR
CXX="${MPICXX}"

AM_PROG_AR
AC_PROG_CPP           dnl sets variable CPP
AC_PROG_RANLIB        dnl sets variable RANLIB
AC_PROG_INSTALL       dnl sets variables INSTALL, INSTALL_DATA, INSTALL_PROGRAM
AC_PROG_LN_S          dnl sets variable LN_S
AC_PROG_MAKE_SET      dnl set $(MAKE) if needed
AC_PROG_RANLIB        dnl set variable RANLIB
AC_CHECK_TOOL([AR],[ar]) dnl set variable AR to appropriate `ar' program

dnl -----------------------------------------------------------------------
AC_OPENMP
AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wpedantic  -Wextra])
AX_APPEND_COMPILE_FLAGS([-pthread])
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}"
LDFLAGS="${LDFLAGS} ${OPENMP_CXXFLAGS}"

dnl -----------------------------------------------------------------------
AX_BLAS([], [AC_MSG_FAILURE([Could not find BLAS library.])])
AX_LAPACK([], [AC_MSG_FAILURE([Could not find LAPACK library.])])

dnl -----------------------------------------------------------------------
dnl CXX checks
dnl
AX_CXX_STATIC_CAST
AX_CXX_TEMPLATES
AX_CXX_TYPENAME
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])
AC_TYPE_SIZE_T

HAVE_SHARED_PTR=no
AC_MSG_CHECKING([for std::shared_ptr])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
      [[#include <memory>
      ]],
      [[
        std::shared_ptr<int> a(new int);
        std::shared_ptr<int> b = a;
      ]]
     )
    ],
    [AC_MSG_RESULT([yes])
     HAVE_SHARED_PTR=std
    ],
    [
      AC_MSG_ERROR([no])
    ]
  )

AC_MSG_CHECKING([for designated initializer])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
      [[#include<cassert>
      ]],
      [[
        struct A { int a, b; };
        A a = {.a = 5, .b = 6};
        assert( a.a == 5 && a.b == 6 );
      ]]
     )
    ],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_ERROR([no])]
  )

dnl -----------------------------------------------------------------------
dnl CHECK CTF
dnl -----------------------------------------------------------------------
if test xYES = x${BUILD_CTF}; then
  AC_MSG_WARN([
               In order to build ctf, type make ctf in the build directory
               ])
else
  CPPFLAGS="$CPPFLAGS -I${LIBCTF_CPATH}"
  LDFLAGS="$LDFLAGS -L${LIBCTF_LD_LIBRARY_PATH} -lctf"
  AC_CHECK_HEADER([ctf.hpp],
                  [],
                  [AC_MSG_ERROR([no ctf.hpp])])
fi

dnl -----------------------------------------------------------------------
dnl CHECK LIBINT
dnl -----------------------------------------------------------------------
if test xYES = x${BUILD_LIBINT}; then
  AC_MSG_WARN([
               In order to build libint, type make libint in the build directory
               ])
else
  CPPFLAGS="$CPPFLAGS -I${LIBINT_CPATH}"
  LDFLAGS="$LDFLAGS -L${LIBINT_LD_LIBRARY_PATH} -lint2"
  AC_CHECK_HEADER([libint2/libint2.hpp],
                  [],
                  [AC_MSG_ERROR([no libint2/libint2.hpp header found])])
fi

dnl -----------------------------------------------------------------------
dnl CHECK EIGEN
dnl -----------------------------------------------------------------------
if test xYES = x${BUILD_EIGEN}; then
  AC_MSG_WARN([
               In order to build libint, type make libint in the build directory
               ])
else
  CPPFLAGS="$CPPFLAGS -I${EIGEN_CPATH}"
  AC_CHECK_HEADER([Eigen/Core],
                  [],
                  [AC_MSG_ERROR([No Eigen/Core header available])])
fi

dnl -----------------------------------------------------------------------
dnl CHECK YAMLCPP
dnl -----------------------------------------------------------------------
if test xYES = x${BUILD_YAMLCPP}; then
  AC_MSG_WARN([
               In order to build yamlcpp, type make yamlcpp in the build directory
               ])
else
  CPPFLAGS="$CPPFLAGS -I${YAMLCPP_CPATH}"
  LDFLAGS="$LDFLAGS -L${YAMLCPP_LD_LIBRARY_PATH} -lyaml-cpp"
  AC_CHECK_HEADER([yaml-cpp/yaml.h],
                  [],
                  [AC_MSG_ERROR([no yaml-cpp/yaml.h header found])])
fi


dnl CUDA NVIDIA -----------------------------------------------------------

AM_CONDITIONAL([WITH_CUDA], [test x${WITH_CUDA} = xyes])
if test x${WITH_CUDA} = xyes; then
  AC_MSG_RESULT([
                     CUDA SUPPORT IS ENABLED
                     -----------------------
  ])
  AC_CHECK_PROGS([NVCC], [nvcc])
  AC_SUBST([CUDA_LDFLAGS])
  AC_DEFINE([HAVE_CUDA],1,[Wether we are using CUDA])
  # TODO: make sure to find cuda and cudart
  # AC_SEARCH_LIBS([cudaMalloc], [cuda cudart], [FOUND_CUDA=yes])
  CC3S_CUBLAS
fi


dnl -----------------------------------------------------------------------
AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
               ])
AC_CONFIG_HEADERS([config.h])


AC_MSG_RESULT([
                       AUTOMATIC CONFIGURATION
                       -----------------------
                       Vendor   = $ax_cv_cxx_compiler_vendor
                       CXX      = $CXX
                       MPICXX   = $MPICXX
                       CXXFLAGS = $CXXFLAGS
                       CPPFLAGS = $CPPFLAGS
                       LDFLAGS  = $LDFLAGS
                       BLAS = ${BLAS_LIBS}
                       LIBS = ${LIBS}
])

AC_OUTPUT
