# SYNOPSIS  -*- mode: autoconf; -*-
#
# 	SISI_SCALAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
# 	This macro looks for a library that implements the SCALAPACK linear-algebra
# 	interface (see http://www.netlib.org/scalapack/). On success, it sets the
# 	SCALAPACK_LIBS output variable to hold the requisite library linkages.
#
# 	To link with SCALAPACK, you should link with:
#
# 		$SCALAPACK_LIBS $LIBS
#
# 	in that order. FLIBS is the output variable of the
# 	AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_SCALAPACK), and is
# 	sometimes necessary in order to link with F77 libraries. Users will also
# 	need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
# 	reason.
#
# 	Many libraries are searched for, from ATLAS to CXML to ESSL. The user
# 	may also use --with-scalapack=<lib> in order to use some specific SCALAPACK
# 	library <lib>. In order to link successfully, however, be aware that you
# 	will probably need to use the same Fortran compiler (which can be set
# 	via the F77 env. var.) as was used to compile the SCALAPACK library.
#
# 	ACTION-IF-FOUND is a list of shell commands to run if a SCALAPACK library is
# 	found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
# 	not found. If ACTION-IF-FOUND is not specified, the default action will
# 	define HAVE_SCALAPACK.
#
# LICENSE
#
# 	Copyright (c) 2022 Alejandro Gallo <aamsgallo@gmail.com>
#
# 	This program is free software: you can redistribute it and/or modify it
# 	under the terms of the GNU General Public License as published by the
# 	Free Software Foundation, either version 3 of the License, or (at your
# 	option) any later version.
#
# 	This program is distributed in the hope that it will be useful, but
# 	WITHOUT ANY WARRANTY; without even the implied warranty of
# 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# 	Public License for more details.
#
# 	You should have received a copy of the GNU General Public License along
# 	with this program. If not, see <https://www.gnu.org/licenses/>.
#
# 	As a special exception, the respective Autoconf Macro's copyright owner
# 	gives unlimited permission to copy, distribute and modify the configure
# 	scripts that are the output of Autoconf when processing the Macro. You
# 	need not follow the terms of the GNU General Public License when using
# 	or distributing such scripts, even though portions of the text of the
# 	Macro appear in them. The GNU General Public License (GPL) does govern
# 	all other use of the material that constitutes the Autoconf Macro.
#
# 	This special exception to the GPL applies to versions of the Autoconf
# 	Macro released by the Autoconf Archive. When you make and distribute a
# 	modified version of the Autoconf Macro, you may extend this special
# 	exception to the GPL to apply to your modified version as well.

AC_DEFUN([SISI_SCALAPACK], [
AC_PREREQ([2.55])
AC_REQUIRE([AC_CANONICAL_HOST])
AC_REQUIRE([AX_BLAS])

AC_ARG_WITH(scalapack,
            [AS_HELP_STRING([--with-scalapack=<lib>],
                            [use SCALAPACK library <lib>|auto|no])])

ax_scalapack_ok=no

case $with_scalapack in
    auto | "")
        ax_scalapack_ok=auto
    ;;
    no)
        ax_scalapack_ok=disable
        ;;
    -* | */* | *.a | *.so | *.so.* | *.dylib | *.dylib.* | *.o)
        SCALAPACK_LIBS="$with_scalapack"
        ;;
    *)
        SCALAPACK_LIBS="-l$with_scalapack"
        ;;
esac

dnl Get fortran linker names of SCALAPACK functions to check for.
AC_F77_FUNC(pcheev)

mkl_scalapack_lib="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"

## Common suffixes for better discoverability
sisi_scalapack_common_suffixes="-lgfortran"

ax_blas_save_LIBS="$LIBS"

################################################################################
# $SCALAPACK_LIBS testing
#
## First, check SCALAPACK_LIBS environment variable
## And go through the possible preffixes
save_SCALAPACK_LIBS="$SCALAPACK_LIBS"
for suffix in "" $sisi_scalapack_common_suffixes; do

  if test "x$save_SCALAPACK_LIBS" != x; then
      save_LIBS="$LIBS"
      test_SCALAPACK_LIBS="$save_SCALAPACK_LIBS $suffix"
      LIBS="$test_SCALAPACK_LIBS $LIBS"
      AC_MSG_CHECKING([for $pcheev in $test_SCALAPACK_LIBS])
      AC_TRY_LINK_FUNC([$pcheev],
                       [ax_scalapack_ok=yes
                        SCALAPACK_LIBS="$test_SCALAPACK_LIBS"],
                       [SCALAPACK_LIBS=""
                        ax_scalapack_ok=no])
      AC_MSG_RESULT([$ax_scalapack_ok])
      LIBS="$save_LIBS"
  fi

done

################################################################################
# MKL TESTING
#
if test "x$ax_scalapack_ok" != xyes ; then
if test -d "${MKLROOT}" ; then
    save_LIBS="$LIBS"
    LIBS="$mkl_scalapack_lib $LIBS"
    AC_MSG_CHECKING([for $pcheev in $mkl_scalapack_lib])
    AC_TRY_LINK_FUNC([$pcheev],
                     [ax_scalapack_ok=yes
                      SCALAPACK_LIBS="$mkl_scalapack_lib"],
                     [SCALAPACK_LIBS=""])
    AC_MSG_RESULT($ax_scalapack_ok)
    LIBS="$save_LIBS"
fi
fi

################################################################################
# -LSCALAPACK TESTING
#
for suffix in "" $sisi_scalapack_common_suffixes; do
  if test "x$ax_scalapack_ok" != xyes ; then

      save_LIBS="$LIBS"
      SCALAPACK_LIBS="-lscalapack $suffix"
      LIBS="$SCALAPACK_LIBS $LIBS"
      AC_MSG_CHECKING([for $pcheev in $SCALAPACK_LIBS])
      AC_TRY_LINK_FUNC([$pcheev],
                       [ax_scalapack_ok=yes],
                       [SCALAPACK_LIBS=""
                        ax_scalapack_ok=no])
      AC_MSG_RESULT([$ax_scalapack_ok])
      LIBS="$save_LIBS"

  fi
done

################################################################################
# DEFAULT TESTING
#
for suffix in "" $sisi_scalapack_common_suffixes; do
if test "x$ax_scalapack_ok" != xyes ; then
    save_LIBS="$LIBS"
    LIBS="$suffix $LIBS"
    AC_MSG_CHECKING([for $pcheev in $LIBS])
    AC_TRY_LINK_FUNC([$pcheev],
                     [ax_scalapack_ok=yes
                      SCALAPACK_LIBS="$suffix"],
                     [SCALAPACK_LIBS=""
                      ax_scalapack_ok=no])
    AC_MSG_RESULT([$ax_scalapack_ok])
    LIBS="$save_LIBS"
fi
done


if test x"$ax_scalapack_ok" == xyes; then
    ifelse([$1],
           ,
           AC_DEFINE([HAVE_SCALAPACK],
                     [1],
                     [Define if you have a SCALAPACK library.]),
           [$1])
    :
else
    ax_scalapack_ok=no
    $2
fi

])dnl DEFUN
