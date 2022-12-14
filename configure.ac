dnl  This is an autoconf script.
dnl  To rebuild the `configure' script from this, execute the command
dnl 	autoconf
dnl  in the directory containing this script.
dnl

AC_INIT([sisi4s],[0.0.0],[https://gihtub.com/alejandrogallo/sisi4s])
AC_CONFIG_SRCDIR([configure.ac])
AC_PREREQ(2.69)
AC_CONFIG_MACRO_DIRS([./etc/m4])
AC_CONFIG_AUX_DIR([./etc/autotools])
AM_INIT_AUTOMAKE([foreign -Wall])

dnl -----------------------------------------------------------------------
dnl External dependencies
dnl -----------------------------------------------------------------------
EXTERNAL_DEPENDENCIES=""
EXTERNAL_LDFLAGS=""
EXTERNAL_CPPFLAGS=""
EXTERNAL_PATH="\${abs_top_builddir}"

external_ldflags_add() {
    EXTERNAL_LDFLAGS="${2}${EXTERNAL_PATH}/extern/$1 $3 $EXTERNAL_LDFLAGS"
}

external_cppflags_add() {
    EXTERNAL_CPPFLAGS="-I${EXTERNAL_PATH}/extern/$1 $EXTERNAL_CPPFLAGS"
}

AC_SUBST([EXTERNAL_DEPENDENCIES])
AC_SUBST([EXTERNAL_LDFLAGS])
AC_SUBST([EXTERNAL_CPPFLAGS])
AC_SUBST([EXTERNAL_PATH])

EXTERN_BUILD_PATH="${ac_pwd}/extern/build/"
EXTERN_SRC_PATH="${ac_pwd}/extern/src/"
AC_SUBST([EXTERN_BUILD_PATH])
AC_SUBST([EXTERN_SRC_PATH])

CTF_COMMIT="53ae5daad851bf3b198ebe1fa761c13b12291116"
YAML_COMMIT="c9460110e072df84b7dee3eb651f2ec5df75fb18"
LIBINT_COMMIT="v2.6.0"
EIGEN_BRANCH="3.3.7"
AC_SUBST([CTF_COMMIT])
AC_SUBST([YAML_COMMIT])
AC_SUBST([LIBINT_COMMIT])
AC_SUBST([EIGEN_BRANCH])

dnl -----------------------------------------------------------------------
dnl Configure switches

dnl CTF WITH OPTION

AC_ARG_WITH(ctf,
            [AS_HELP_STRING([--with-ctf],
                           [prefix for CTF includes and libraries] )],
            [LIBCTF_PATH="`readlink -f $withval`";
             LIBCTF_CPATH="`readlink -f $withval`/include";
             LIBCTF_LD_LIBRARY_PATH="`readlink -f $withval`/lib";
             LIBS="${LIBS} -L${LIBCTF_LD_LIBRARY_PATH} -lctf"
             BUILD_CTF=NO],
            [BUILD_CTF=YES
             external_ldflags_add "ctf/lib/libctf.a"
             external_cppflags_add "ctf/include/"
             EXTERNAL_DEPENDENCIES="ctf ${EXTERNAL_DEPENDENCIES}"])

dnl LIBINT WITH OPTION
AC_ARG_WITH(libint,
            [AS_HELP_STRING([--with-libint],
                            [prefix for libint includes and libraries] )],
            [LIBINT_PATH="`readlink -f $withval`";
             LIBINT_CPATH="`readlink -f $withval`/include";
             LIBINT_LD_LIBRARY_PATH="`readlink -f $withval`/lib";
             LIBS="${LIBS} -L${LIBINT_LD_LIBRARY_PATH} -lint2"
             BUILD_LIBINT=NO],
            [BUILD_LIBINT=YES])

AC_ARG_ENABLE([libint],
              [AS_HELP_STRING([--disable-libint],
                              [Don't enable libint if it is not going to be used])],
              [DISABLE_LIBINT=yes
               BUILD_LIBINT=NO],
              [DISABLE_LIBINT=no
               AC_DEFINE([HAVE_LIBINT],
                         [1],
                         [Wether libint is used])])

# Do this here after --disable-libint has been processed
if test x$BUILD_LIBINT = xYES; then
  EXTERNAL_DEPENDENCIES="libint ${EXTERNAL_DEPENDENCIES}"
  external_ldflags_add "libint/lib/libint2.a"
  external_cppflags_add "libint/include/"
fi

dnl EIGEN WITH OPTION
AC_ARG_WITH(eigen,
            [AS_HELP_STRING([--with-eigen],
                            [prefix for libint includes and libraries] )],
            [EIGEN_PATH="`readlink -f $withval`";
             EIGEN_CPATH="`readlink -f $withval`/include/eigen3";
             BUILD_EIGEN=NO],
            [BUILD_EIGEN=YES
             external_cppflags_add "eigen/include/eigen3/"
             EXTERNAL_DEPENDENCIES="eigen ${EXTERNAL_DEPENDENCIES}"])

dnl LIBYAML WITH OPTION
AC_ARG_WITH(yamlcpp,
            [AS_HELP_STRING([--with-yamlcpp],
                            [prefix for yamlcpp includes and libraries] )],
            [YAMLCPP_PATH="`readlink -f $withval`";
             YAMLCPP_CPATH="`readlink -f $withval`/include";
             YAMLCPP_LD_LIBRARY_PATH="`readlink -f $withval`/lib";
             BUILD_YAML=NO],
            [BUILD_YAML=YES
             external_cppflags_add "yaml-cpp/include/"
             external_ldflags_add "yaml-cpp/lib/libyaml-cpp.a"
             EXTERNAL_DEPENDENCIES="yaml-cpp ${EXTERNAL_DEPENDENCIES}"])

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

AC_ARG_ENABLE([docs],
              [AS_HELP_STRING([--enable-docs],
                              [Enable building documentation])],
                              [build_docs=yes],
                              [build_docs=no])

AC_ARG_VAR([NVCC],
           [Path to the nvidia cuda compiler.])
AC_ARG_VAR([CUDA_LDFLAGS],
           [LDFLAGS to find libraries -lcuda, -lcudart, -lcublas.])
AC_ARG_VAR([CUDA_CXXFLAGS],
           [CXXFLAGS to find the CUDA headers])


dnl -----------------------------------------------------------------------

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

CXX_VENDOR=$ax_cv_cxx_compiler_vendor
CXX="${MPICXX}"

AM_PROG_AR
AC_PROG_CPP              dnl sets variable CPP
AC_PROG_RANLIB           dnl sets variable RANLIB
AC_PROG_INSTALL          dnl sets variables INSTALL, INSTALL_DATA,
                         dnl                INSTALL_PROGRAM
AC_PROG_LN_S             dnl sets variable LN_S
AC_PROG_AWK
AC_PROG_GREP
AC_PROG_INSTALL
AC_PROG_SED
AC_PROG_MAKE_SET         dnl set $(MAKE) if needed
AC_PROG_RANLIB           dnl set variable RANLIB
AC_CHECK_TOOL([TR],[tr])

dnl -----------------------------------------------------------------------
AC_OPENMP
AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wpedantic  -Wextra])
AX_APPEND_COMPILE_FLAGS([-pthread])
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}"
LDFLAGS="${LDFLAGS} ${OPENMP_CXXFLAGS}"

dnl -----------------------------------------------------------------------

AC_MSG_RESULT([
                    LINEAR ALGEBRA LIBRARIES
                    ------------------------
])

AX_BLAS([LIBS="${LIBS} ${BLAS_LIBS}"],
        [AC_MSG_FAILURE([Could not find BLAS library.])])
AX_LAPACK([LIBS="${LIBS} ${LAPACK_LIBS}"],
          [AC_MSG_FAILURE([Could not find LAPACK library.])])

SISI_SCALAPACK([LIBS="${LIBS} ${SCALAPACK_LIBS}"],
               [AC_MSG_FAILURE([could not find a working scalapack installation.])])

dnl -----------------------------------------------------------------------
dnl CXX checks
dnl
AC_MSG_RESULT([
                          C++ features test
                          -----------------
])
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
AM_CONDITIONAL([WITH_BUILD_CTF], [test x${BUILD_CTF} = xYES])
if test xNO = x${BUILD_CTF}; then
  CPPFLAGS="$CPPFLAGS -I${LIBCTF_CPATH}"
  AC_CHECK_HEADER([ctf.hpp],
                  [],
                  [AC_MSG_ERROR([no ctf.hpp])])
fi

dnl -----------------------------------------------------------------------
dnl CHECK LIBINT
dnl -----------------------------------------------------------------------
AM_CONDITIONAL([WITH_BUILD_LIBINT], [test x${BUILD_LIBINT} = xYES])
AM_CONDITIONAL([DISABLE_LIBINT], [test x${DISABLE_LIBINT} = xyes])
LIBINT_MAX_AM=4
AC_ARG_VAR([LIBINT_MAX_AM],
           [Maximum angular momentum for libint. Default: 4],)
if test xNO = x${BUILD_LIBINT} && test x${DISABLE_LIBINT} = xno; then
  CPPFLAGS="$CPPFLAGS -I${LIBINT_CPATH}"
  #LDFLAGS="$LDFLAGS -L${LIBINT_LD_LIBRARY_PATH} -lint2"
  AC_CHECK_HEADER([libint2.h],
                  [],
                  [AC_MSG_ERROR([no libint2.hpp header found in $CPPFLAGS])])
fi


dnl -----------------------------------------------------------------------
dnl CHECK EIGEN
dnl -----------------------------------------------------------------------
AM_CONDITIONAL([WITH_BUILD_EIGEN], [test x${BUILD_EIGEN} = xYES])
if test xNO = x${BUILD_EIGEN}; then
  CPPFLAGS="$CPPFLAGS -I${EIGEN_CPATH}"
  AC_CHECK_HEADER([Eigen/Core],
                  [],
                  [AC_MSG_ERROR([No Eigen/Core header available with $EIGEN_CPATH])])
fi

dnl -----------------------------------------------------------------------
dnl CHECK YAMLCPP
dnl -----------------------------------------------------------------------
AM_CONDITIONAL([WITH_BUILD_YAML], [test x${BUILD_YAML} = xYES])
if test xNO = x${BUILD_YAML}; then
  CPPFLAGS="$CPPFLAGS -I${YAMLCPP_CPATH}"
  #LDFLAGS="$LDFLAGS -L${YAMLCPP_LD_LIBRARY_PATH} -lyaml-cpp"
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
dnl Docs
AM_CONDITIONAL([WITH_DOCS], [test x${build_docs} = xyes])
if test x${build_docs} == xyes; then
  AC_CHECK_PROGS([DOXYGEN], [doxygen])
  AC_CHECK_PROGS([DOT], [dot])
  AX_PROG_EMACS
  AX_ELISP_CHECK_FEATURE([org])
  AX_ELISP_CHECK_FBOUNDP([package-initialize], [package])
  AX_ELISP_CHECK_FBOUNDP([package-install], [package])
  AX_ELISP_CHECK_FBOUNDP([string-join], [subr-x])
  AX_ELISP_CHECK_BOUNDP([org-babel-post-tangle-hook], [ob])
fi

AC_DEFINE_UNQUOTED([SISI_COMMIT],
                   ["$(git describe --always --dirty)"],
                   [The current git status at compile time])

AC_DEFINE_UNQUOTED([COMPILER_VERSION],
                   ["$($MPICXX --version | $TR "\n" " ")"],
                   [Compiler version])

dnl -----------------------------------------------------------------------
AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
                 docs/Makefile
                 docs/Doxyfile
                 docs/config.el
                 docs/conf.py
                 integration-tests/config.mk
                 integration-tests/Makefile
                 extern/Makefile
                 extern/yaml-cpp/Makefile
                 extern/eigen/Makefile
                 extern/libint/Makefile
                 extern/ctf/Makefile
                 .dir-locals.el
               ])
AC_CONFIG_HEADERS([config.h])


AC_MSG_RESULT([
                       AUTOMATIC CONFIGURATION
                       -----------------------
                       Vendor  =  $ax_cv_cxx_compiler_vendor
                       CXX=$CXX
                       MPICXX=$MPICXX
                       CXXFLAGS=$CXXFLAGS
                       CPPFLAGS=$CPPFLAGS
                       LDFLAGS=$LDFLAGS
                       BLAS_LIBS=${BLAS_LIBS}
                       SCALAPACK_LIBS=${SCALAPACK_LIBS}
                       LIBS=${LIBS}
])

if test -n "${EXTERNAL_DEPENDENCIES}"; then
  AC_MSG_RESULT([
                 In order to build the project you need first to
                 make the dependencies by doing

                   make extern

                 or alternatively
                 ])
  for extern in ${EXTERNAL_DEPENDENCIES}; do
  AC_MSG_RESULT([                 - make -C extern/$extern])
  done
fi

AC_OUTPUT
