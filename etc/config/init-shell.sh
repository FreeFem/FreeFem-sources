ffapp=@ff_prefix@
version=@VERSION@
ff_prefix_dir_lib=@ff_prefix_dir_lib@
ff_prefix_dir_lib_mpi=@ff_prefix_dir_lib_mpi@
ffpp=`which FreeFem++`
ffsrc=@abs_top_srcdir@
fftmp=/tmp/ff-$$
bindir=@bindir@
if [ "$ffpp" != "$ffapp/bin/FreeFem++" ] ; then
  echo sorry appli must be install !!
  echo $ffpp 
  echo $ffapp/bin/FreeFem++
  exit 1
fi
arch=xxx
case $machine in
arm64e)
 arch=M1;;
x86_64)
  arch=intel;;
*)
  echo unknow machine $machine
  exit 1;;
esac
echo " arch = $arch"

