// Include a platform-specific configuration file
// ----------------------------------------------

// $Id$

// This wrapper is necessary for platforms where the configure script
// does not run.

#ifndef CONFIG_WRAPPER_H
#define CONFIG_WRAPPER_H

#ifdef __MWERKS__
#include <config-macos9.h>
#else
#include <config.h>
#endif

#endif // CONFIG_WRAPPER_H
