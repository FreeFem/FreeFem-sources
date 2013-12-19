*** /home/alh/ff/draft/download/blas/Makefile.am	2013-12-19 11:25:49.344446251 +0100
--- /tmp/Makefile.am~other.wWaife	2013-12-19 11:25:49.344446251 +0100
***************
*** 25,41 ****
  BLAS_TGZ=$(DIRPKG)/blas.tgz
  CBLAS_TGZ=$(DIRPKG)/cblas.tgz
  
! COMMON_PACKTITLE=BLAS
! COMMON_PACKAGES=$(BLAS_TGZ) $(CBLAS_TGZ) ../pkg/OpenBLAS.tar.gz
! include ../common.mak
  
  # Downloading and compiling the Generic Blas
  # ------------------------------------------
  
  noinst_LIBRARIES=@DOWNLOADED_BLAS@
- compile::$(noinst_LIBRARIES)
  EXTRA_LIBRARIES=libf77blas.a libcblas.a
- compile::@DOWNLOADED_BLAS@
  
  # List of files to compile (do not list them in *_SOURCES to prevent them from being included in distributions).
  
--- 25,38 ----
  BLAS_TGZ=$(DIRPKG)/blas.tgz
  CBLAS_TGZ=$(DIRPKG)/cblas.tgz
  
! EXTRA_DIST= \
! openblas.patches 
  
  # Downloading and compiling the Generic Blas
  # ------------------------------------------
  
  noinst_LIBRARIES=@DOWNLOADED_BLAS@
  EXTRA_LIBRARIES=libf77blas.a libcblas.a
  
  # List of files to compile (do not list them in *_SOURCES to prevent them from being included in distributions).
  
***************
*** 122,137 ****
  	touch CBLAS/fait
  $(CBLAS_SOURCES): CBLAS/fait
  
  clean-local::
  	-rm -r BLAS CBLAS blas_xerbla.f
  
! # <<OpenBLAS>> ALH - 18/9/13 - Downloading and building the OpenBLAS
! # ------------------------------------------------------------------
  
  # to activate this, see [[file:../../configure.ac::OpenBLAS]]
  
  all-local::@COMPILE_OPENBLAS@
- compile::@COMPILE_OPENBLAS@
  generic:
  openblas:links.done
  
--- 119,142 ----
  	touch CBLAS/fait
  $(CBLAS_SOURCES): CBLAS/fait
  
+ $(BLAS_TGZ):
+ 	mkdir -p $(DIRPKG)
+ 	cd $(DIRPKG) && @WGET@ -N http://www.netlib.org/blas/blas.tgz
+ 
+ $(CBLAS_TGZ):
+ 	mkdir -p $(DIRPKG)
+ 	cd $(DIRPKG) && @WGET@ -N http://www.netlib.org/blas/blast-forum/cblas.tgz
+ 
+ 
  clean-local::
  	-rm -r BLAS CBLAS blas_xerbla.f
  
! # ALH - 18/9/13 - Downloading and building the OpenBLAS
! # -----------------------------------------------------
  
  # to activate this, see [[file:../../configure.ac::OpenBLAS]]
  
  all-local::@COMPILE_OPENBLAS@
  generic:
  openblas:links.done
  
***************
*** 158,178 ****
  
  # The OpenBLAS directory is updated during the compilation, so the patching step should not depend on the directory
  # date
- 
  openpatches.done:opendownload.done
- 	tar xvzf ../pkg/OpenBLAS.tar.gz
  	patch -u -p1 < openblas.patches
  	touch $@
- 
- download::opendownload.done
  opendownload.done:
! 	../getall -o OpenBLAS
  	touch $@
  
! # ALH - 6/11/13 - Since OpenBLAS is quite long to compile and it does not change very often, it's only cleaned as part
! # of the specific target 'veryclean'.
  
! veryclean::
  	-rm -r *.done OpenBLAS
  
  # Local Variables:
--- 163,179 ----
  
  # The OpenBLAS directory is updated during the compilation, so the patching step should not depend on the directory
  # date
  openpatches.done:opendownload.done
  	patch -u -p1 < openblas.patches
  	touch $@
  opendownload.done:
! 	git clone https://github.com/xianyi/OpenBLAS.git
  	touch $@
  
! # Since some libraries (VTK, FLTK, Openblas) are quite long to compile and they do not change very often (we don't
! # usually change them in FFCS), they are only cleaned with specific target 'allclean'.
  
! allclean::
  	-rm -r *.done OpenBLAS
  
  # Local Variables:
