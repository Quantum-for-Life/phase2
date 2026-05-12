#ifndef PHASE2_DET_SMALL_H
#define PHASE2_DET_SMALL_H

#include <stdint.h>

/*
 * det_small - small dense determinant via LU with partial pivoting.
 *
 * Public API:
 *   double det_small(const double *A, uint32_t n);
 *
 * Operates on a row-major n x n real matrix.  Returns the determinant.
 * No external libraries; heap-free; bounded for n <= DET_SMALL_MAX_N.
 * For larger n, returns 0.0 (caller must check the bound).
 *
 * The routine is intended for the Slater-Condon expansion
 * (`state_prep_coeff_expand`), where n equals the alpha/beta
 * occupation count (n <= 17 at half-filling N = 34, tapered).
 */

#define DET_SMALL_MAX_N (32u)

double det_small(const double *A, uint32_t n);

#endif /* PHASE2_DET_SMALL_H */
