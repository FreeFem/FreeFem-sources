#ifdef __GNUC__ 
//#if (__GNUC_MINOR__ < 90) && (D__GNUG__ <= 3) 
#define INCLUDE_TEMPLATE_DEFINITION
//#endif
#endif
#ifdef XLC_TEMPINC
#define INCLUDE_TEMPLATE_DEFINITION
#endif
#ifdef __MWERKS__
#define INCLUDE_TEMPLATE_DEFINITION
#endif
#ifdef __BCPLUSPLUS__
#define INCLUDE_TEMPLATE_DEFINITION
#endif
#ifdef _MSC_VER
#define INCLUDE_TEMPLATE_DEFINITION
#endif
