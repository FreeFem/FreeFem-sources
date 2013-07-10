Summary: FreeFem++
Name: freefem++
Version: 2.8
Release: 0
Source: %{name}-%{version}.tar.gz
Patch:   %{name}-config.patch
Patch1:   %{name}-gcc4.patch
%if %{?_with_cadna:1}%{!?_with_cadna:0} 
Source2: CadnaC_gcc-3.2_Linux_i386.tar.gz
Patch2:  cadna-gcc4.patch
%endif
License: GPL
Group: Applications/Engineering
URL: http://www.freefem.org/ff++/
Packager: Christophe  Trophime <christophe.trophime@grenoble.cnrs.fr>
Prereq: /sbin/install-info
Buildroot: %{_tmppath}/%{name}-buildroot
Requires: arpack, ufsparse
BuildRequires: arpack-devel, ufsparse-devel
BuildRequires: fltk >= 1.1
BuildRequires: fltk >= 1.1
BuildRequires: gsl >= 1.2
BuildRequires: rpm >= 4.1
%if %{?_with_mpi:1}%{!?_with_mpi:0} 
BuildRequires: lam
%endif
%{!?_without_freedesktop:BuildRequires: desktop-file-utils}
Requires: mesa-libGL >= 6.7.0-9
Requires: mesa-libGLU >= 6.7.0-9
Requires: gsl >= 1.2
Requires: fltk >= 1.1
Prefix: /usr

%description 
FreeFem++ is an implementation of a language dedicated to the finite element method. 
It enables you to solve Partial Differential Equations (PDE) easily.

Problems involving PDE from several branches of physics such as fluid-structure interactions 
require interpolations of data on several meshes and their manipulation within one program. 
FreeFem++ includes a fast quadtree-based interpolation algorithm and a language for the manipulation 
of data on multiple meshes (generated with bamg).



%prep

%setup -q -n %{name}-%{version}
%patch -p1 -b .umfpack
%patch1 -p1 -b .gcc4

%if %{?_with_cadna:1}%{!?_with_cadna:0}
mkdir -p cadna
mkdir -p download/cadna
tar zxvf %{SOURCE2} -C cadna 
mv cadna/include/cadnafree.h download/cadna
mv cadna/lib/libcadnafreeC.a download/cadna
pushd download/cadna/
ln -sf libcadnafreeC.a libcadnafree.a
popd
%patch2 -p1 -b .cadna-gcc4.patch
%endif

autoreconf -f -i

%build
%if %{?_with_mpi:1}%{!?_with_mpi:0}
%configure --with-mpi=lam
%else
%configure --without-mpi --with-blas="-L/usr/lib/atlas -lf77blas -lcblas"
%endif
make

%install
rm -rf $RPM_BUILD_ROOT

make install DESTDIR=$RPM_BUILD_ROOT



%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{_bindir}/FreeFem++
%{_bindir}/FreeFem++-cs
%{_bindir}/FreeFem++-nw
%{_bindir}/FreeFem++-glx
%{_bindir}/FreeFem++-ide
%{_bindir}/FreeFem++-server
%{_bindir}/FreeFem++-client
%{_bindir}/bamg
%{_bindir}/cvmsh2
%{_bindir}/drawbdmesh
%if %{?_with_mpi:1}%{!?_with_mpi:0}
%{_bindir}/FreeFem++-mpi
%endif

