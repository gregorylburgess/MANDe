AC_INIT([ncdf], 1.6.5)

AC_ARG_WITH([nc-config],
    AC_HELP_STRING([--with-nc-configh=PATH], [location of nc-config]),
    [NC_CONFIG=$withval])
if test -z "${NC_CONFIG}"; then
   AC_PATH_PROG(NC_CONFIG, nc-config)
else
  if ! test -x "${NC_CONFIG}"; then
    NC_CONFIG=""
    AC_MSG_WARN(['NC_CONFIG' is set to an invalid value and will be ignored])
  fi
fi

if test -n "${NC_CONFIG}"; then
  PKG_CPPFLAGS=`${NC_CONFIG} --cflags`
  PKG_LIBS=`${NC_CONFIG} --libs`
else
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
        echo "Could not determine R_HOME"
	exit 1
fi
CC=`${R_HOME}/bin/R CMD config CC`
CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
AC_PROG_CC
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CPP

AC_ARG_WITH([netcdf-include],
    AC_HELP_STRING([--with-netcdf-include=DIR], 
                   [location of netcdf header files]),
    [netcdf_include_path=$withval])
if test -n "${netcdf_include_path}"; then
   CPPFLAGS="-I${netcdf_include_path} ${CPPFLAGS}"
elif test -n ${NETCDF_INCLUDE}; then
   CPPFLAGS="-I${NETCDF_INCLUDE} ${CPPFLAGS}"
fi
AC_CHECK_HEADERS(netcdf.h, , AC_MSG_ERROR(netcdf header netcdf.h not found))

AC_ARG_WITH([netcdf-lib],
    AC_HELP_STRING([--with-netcdf-lib=DIR], [location of netcdf libraries]),
    [netcdf_lib_path=$withval])
if test -n "${netcdf_lib_path}"; then
   LIBS="-L${netcdf_lib_path}"
elif test -n "${NETLIB_LIB}"; then
   LIBS="-L${NETCDF_LIB}"
fi

dnl Extra libraries: we try '-lhdf5 -lhdf5_hl -lcurl -lz' if needed
AC_ARG_WITH([netcdf-extra-libs],
    AC_HELP_STRING([--with-netcdf-extra-libs=], [extra libs for netcdf]),
    [netcdf_extras=$withval])
if test -n "${netcdf_extras}"; then
   EXTRAS="-L${netcdf_extras}"
elif test -n "${NETLIB_EXTRAS}"; then
   EXTRAS="-L${NETCDF_EXTRAS}"
fi

# Optional library path for HDF5
AC_ARG_WITH([hdf5-lib],
    AC_HELP_STRING([--with-hdf5-lib=DIR], [location of hdf5 libraries]),
    [hdf5_lib_path=$withval])
if test -n "$hdf5_lib_path"; then
   LIBS=-L${hdf5_lib_path}"
eiif test -n "$HDF5_LIB"; then
   LIBS=-L${HDF5_LIB}"
fi


AC_CHECK_LIB(netcdf, nc_open, [],
[
    AC_CHECK_LIB(netcdf, nc_create, [have_netcdf_hdf5=TRUE],
                 AC_MSG_ERROR(netcdf library not found),
                 [-lhdf5 -lhdf5_hl -lcurl -lz])
    if test "${have_netcdf_hdf5}" = TRUE; then
        LIBS="${LIBS} -lnetcdf -lhdf5 -lhdf5_hl -lcurl -lz"
    fi
], [$EXTRAS])

PKG_CPPFLAGS=${CPPFLAGS}
PKG_LIBS=${LIBS}
fi

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)

AC_OUTPUT(src/Makevars)
