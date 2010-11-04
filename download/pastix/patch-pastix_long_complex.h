--- pastix.h	2010-06-30 11:35:14.000000000 +0200
+++ ../../../pastix/pastix_long_complex.h	2010-06-28 17:53:56.000000000 +0200
@@ -1,9 +1,16 @@
+#ifdef __cplusplus
+#include <complex>
+#define pastix_float_t   std::complex<double> 
+extern "C" {
+  MPI_Datatype GetMpiType() ;
+#else
 #include <complex.h>
+#define pastix_float_t   double complex
+#endif
 
 #define pastix_int_t     long
 #define pastix_uint_t    unsigned long
 #define MPI_PASTIX_INT   MPI_LONG
-#define pastix_float_t   double complex
 #define MPI_PASTIX_FLOAT GetMpiType()
 #define INT              pastix_int_t
 #define UINT             pastix_uint_t
@@ -744,3 +751,7 @@
 #define MTX_ISRHS(a) (a)[0]!='\0'
 
 /* **************************************** */
+#ifdef __cplusplus
+}
+#endif
+
