FILES = \
	ax_blas.m4 \
	ax_cxx_compile_stdcxx.m4 \
	ax_cxx_static_cast.m4 \
	ax_lapack.m4 \
	ax_mpi.m4 \
	ax_cxx_typename.m4 \
	ax_cxx_templates.m4 \
	ax_openmp.m4 \
	ax_prog_cxx_mpi.m4 \
	ax_prog_emacs.m4 \
	ax_elisp.m4 \


update: autoconf-archive
	for f in $(FILES); do cp -v $</m4/$$f $$f; done

autoconf-archive:
	git clone git://git.sv.gnu.org/autoconf-archive.git
