--- pastix_release_2200/install/pastix_int_complex.h	2015-03-11 22:32:55.000000000 +0100
+++ pastix_int_complex.h	2015-03-11 22:38:41.000000000 +0100
@@ -1,9 +1,18 @@
+/*  patched  by  F. Hecht for freefem++  */ 
+#ifdef __cplusplus
+#include <complex>
+#define pastix_float_t   std::complex<double> 
+extern "C" {
+  MPI_Datatype GetMpiType() ;
+#else
 #include <complex.h>
+#define pastix_float_t   double complex
+#endif
+
 
 #define pastix_int_t     int
 #define pastix_uint_t    unsigned int
 #define MPI_PASTIX_INT   MPI_INT
-#define pastix_float_t   double complex
 #define MPI_PASTIX_FLOAT GetMpiType()
 #define INT              pastix_int_t
 #define UINT             pastix_uint_t
@@ -744,3 +753,6 @@
 #define MTX_ISRHS(a) (a)[0]!='\0'
 
 /* **************************************** */
+#ifdef __cplusplus
+}
+#endif
