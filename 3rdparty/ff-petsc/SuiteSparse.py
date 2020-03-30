import config.package

class Configure(config.package.Package):
  def __init__(self, framework):
    config.package.Package.__init__(self,framework)
    self.version           = '5.7.1'
    self.versionname       = 'UMFPACK_MAIN_VERSION.UMFPACK_SUB_VERSION.UMFPACK_SUBSUB_VERSION'
    self.versioninclude    = 'umfpack.h'
    #  Note that there is not SuitSparse version number in the code, the only version information is for UMFPACK
    self.download          = ['http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.3.tar.gz',
                              'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/SuiteSparse-4.4.3.tar.gz']
    self.liblist           = [['libumfpack.a','libklu.a','libcholmod.a','libbtf.a','libccolamd.a','libcolamd.a','libcamd.a','libamd.a','libsuitesparseconfig.a'],
                             ['libumfpack.a','libklu.a','libcholmod.a','libbtf.a','libccolamd.a','libcolamd.a','libcamd.a','libamd.a','libsuitesparseconfig.a','librt.a'],
                             ['libumfpack.a','libklu.a','libcholmod.a','libbtf.a','libccolamd.a','libcolamd.a','libcamd.a','libamd.a','libmetis.a','libsuitesparseconfig.a'],
                             ['libumfpack.a','libklu.a','libcholmod.a','libbtf.a','libccolamd.a','libcolamd.a','libcamd.a','libamd.a','libmetis.a','libsuitesparseconfig.a','librt.a']]
    self.functions        = ['umfpack_dl_wsolve','cholmod_l_solve','klu_l_solve']
    self.includes         = ['umfpack.h','cholmod.h','klu.h']
    self.hastests         = 1
    self.hastestsdatafiles= 1
    self.precisions       = ['double']
    return

  def setupHelp(self, help):
    import nargs
    config.package.Package.setupHelp(self, help)
    help.addArgument('SUITESPARSE', '-download-suitesparse-gpu=<bool>',    nargs.ArgBool(None, 0, 'Install SuiteSparse to use GPUs'))
    
  def setupDependencies(self, framework):
    config.package.Package.setupDependencies(self, framework)
    self.blasLapack = framework.require('config.packages.BlasLapack',self)
    self.mathlib    = framework.require('config.packages.mathlib',self)
    self.deps       = [self.blasLapack,self.mathlib]
    if self.argDB['download-suitesparse-gpu']:
      self.cuda       = framework.require('config.packages.cuda',self)
      self.deps.append(self.cuda)
    return

  def Install(self):
    import os
    self.log.write('SuiteSparseDir = '+self.packageDir+' installDir '+self.installDir+'\n')
    if not self.make.haveGNUMake:
      raise RuntimeError('SuiteSparse buildtools require GNUMake. Use --with-make=gmake or --download-make')

    mkfile = 'SuiteSparse_config/SuiteSparse_config.mk'
    g = open(os.path.join(self.packageDir, mkfile), 'w')
    self.setCompilers.pushLanguage('C')
    g.write('CC           = '+self.setCompilers.getCompiler()+'\n')
    g.write('CF           = '+self.removeWarningFlags(self.setCompilers.getCompilerFlags())+'\n')
    self.setCompilers.popLanguage()
    g.write('MAKE         ='+self.make.make+'\n')
    g.write('RANLIB       = '+self.setCompilers.RANLIB+'\n')
    g.write('ARCHIVE      = '+self.setCompilers.AR+' '+self.setCompilers.AR_FLAGS+'\n')
    g.write('RM           = '+self.programs.RM+'\n')
    g.write('MV           = '+self.programs.mv+'\n')
    g.write('CP           = '+self.programs.cp+'\n')
    g.write('CLEAN             = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d\n')
    g.write('INSTALL_LIB       = ' + self.libDir + '\n')
    g.write('INSTALL_INCLUDE   = ' + self.includeDir + '\n')
    if self.blasLapack.mangling == 'underscore':
      flg = ''
    elif self.blasLapack.mangling == 'caps':
      flg = '-DBLAS_CAPS_DOES_NOT_WORK'
    else:
      flg = '-DBLAS_NO_UNDERSCORE'
    g.write('UMFPACK_CONFIG    = '+flg+'\n')
    if self.argDB['download-suitesparse-gpu']:
      if self.defaultIndexSize == 32:
        raise RuntimeError('SuiteSparse only uses GPUs with --with-64-bit-indices')
      if not hasattr(self.compilers, 'CUDAC'):
        raise RuntimeError('Run with --with-cuda to use allow SuiteSparse to compile using CUDA')
      # code taken from cuda.py
      self.pushLanguage('CUDA')
      petscNvcc = self.getCompiler()
      self.popLanguage()
      self.getExecutable(petscNvcc,getFullPath=1,resultName='systemNvcc')
      if hasattr(self,'systemNvcc'):
        nvccDir = os.path.dirname(self.systemNvcc)
        cudaDir = os.path.split(nvccDir)[0]
      else:
        raise RuntimeError('Unable to locate CUDA NVCC compiler')
      g.write('CUDA_ROOT     = '+cudaDir+'\n')
      g.write('GPU_BLAS_PATH = $(CUDA_ROOT)\n')
      g.write('GPU_CONFIG    = -I$(CUDA_ROOT)/include -DGPU_BLAS\n')
# GPU_CONFIG    = -I$(CUDA_ROOT)/include -DGPU_BLAS -DCHOLMOD_OMP_NUM_THREADS=10
      g.write('CUDA_PATH     = $(CUDA_ROOT)\n')
      g.write('CUDART_LIB    = $(CUDA_ROOT)/lib64/libcudart.so\n')
      g.write('CUBLAS_LIB    = $(CUDA_ROOT)/lib64/libcublas.so\n')
      g.write('CUDA_INC_PATH = $(CUDA_ROOT)/include/\n')
      g.write('NV20          = -arch=sm_20 -Xcompiler -fPIC\n')
      g.write('NV30          = -arch=sm_30 -Xcompiler -fPIC\n')
      g.write('NV35          = -arch=sm_35 -Xcompiler -fPIC\n')
      g.write('NVCC          = $(CUDA_ROOT)/bin/nvcc\n')
      g.write('NVCCFLAGS     = -O3 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -Xcompiler -fPIC\n')
      g.write('CHOLMOD_CONFIG    = '+flg+' -DNPARTITION $(GPU_CONFIG)\n')
      self.addDefine('USE_SUITESPARSE_GPU',1)
    else:
      g.write('CHOLMOD_CONFIG    = '+flg+' -DNPARTITION\n')
    g.close()

    if self.installNeeded(mkfile):
      try:
        self.logPrintBox('Compiling and installing SuiteSparse; this may take several minutes')
        self.installDirProvider.printSudoPasswordMessage()
        # SuiteSparse install does not create missing directories, hence we need to create them first 
        output,err,ret = config.package.Package.executeShellCommand(self.installSudo+'mkdir -p '+os.path.join(self.installDir,'lib'), timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand(self.installSudo+'mkdir -p '+os.path.join(self.installDir,'include'), timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/SuiteSparse_config && '+self.make.make+' && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/AMD && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/COLAMD && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/BTF && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/CAMD && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/CCOLAMD && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/CHOLMOD && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/UMFPACK && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)
        output,err,ret = config.package.Package.executeShellCommand('cd '+self.packageDir+'/KLU && '+self.make.make+' library && '+self.installSudo+self.make.make+' install && '+self.make.make+' clean', timeout=2500, log=self.log)

        self.addDefine('HAVE_SUITESPARSE',1)
      except RuntimeError as e:
        raise RuntimeError('Error running make on SuiteSparse: '+str(e))
      self.postInstall(output+err, mkfile)
    return self.installDir

  def consistencyChecks(self):
    config.package.Package.consistencyChecks(self)
    if self.framework.argDB['with-'+self.package] and self.defaultIndexSize == 64 and self.types.sizes['void-p'] == 4:
      raise RuntimeError('SuiteSparse does not support 64bit indices in 32bit (pointer) mode.')
    return
