/*
 * include/fpu.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.  The algorithms in the double-double and
 * quad-double package does not function with the extended mode found in
 * these FPU.
 */
#ifndef _QD_FPU_H
#define _QD_FPU_H

#include <qd/qd_config.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Set the round-to-double flag, and save the old control word in old_cw.
 * If old_cw is NULL, the old control word is not saved.
 */
QD_API void fpu_fix_start(unsigned int *old_cw);

/*
 * Restore the control word.
 */
QD_API void fpu_fix_end(unsigned int *old_cw);

#ifdef __cplusplus
}
#endif

#endif  /* _QD_FPU_H */
