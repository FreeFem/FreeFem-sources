#!C:\msys64\usr\bin\bash.exe --login
source shell mingw64

echo "Job PETsc"

freefemPETScVersion=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile | cut -c9-15)

PETScDir='/c/builds/shared/'

major=$(grep "#define PETSC_VERSION_MAJOR" "$PETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-35 )
minor=$(grep "#define PETSC_VERSION_MINOR" "$PETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-36 )
subminor=$(grep "#define PETSC_VERSION_SUBMINOR" "$PETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-35 )
PETScVersion=$major.$minor.$subminor

# Check version
if [ "$freefemPETScVersion" = "$PETScVersion" ]
then
    echo "PETSc is up-to-date"
    echo "Nothing to do"
else
    echo "PETSc is outdated"
    echo "Update"
    rm -rf $PETScDir/ff-petsc
    
    # configuration & build
    autoreconf -i
    ./configure --enable-generic --enable-optim --enable-download \
        --prefix="$PETScDir"
    ./3rdparty/getall -a -o PETSc
    cd 3rdparty/ff-petsc && make petsc-slepc
fi
