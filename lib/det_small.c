/*
 * det_small - small dense determinant via LU with partial pivoting.
 *
 * Public API:
 *   double det_small(const double *A, uint32_t n);
 *
 * The matrix size is bounded by DET_SMALL_MAX_N; for larger n the
 * routine returns 0.0 and the caller is expected to check the
 * bound.  The implementation copies the input into a fixed
 * stack buffer, performs Gaussian elimination with partial
 * pivoting, and accumulates the product of diagonal pivots
 * times the row-swap sign.
 *
 * Edge cases:
 *   n == 0          returns 1.0  (det of the empty matrix).
 *   any pivot == 0  returns 0.0  (singular).
 */

#include "c23_compat.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "det_small.h"

double det_small(const double *A, const uint32_t n)
{
	if (n > DET_SMALL_MAX_N)
		return 0.0;
	if (n == 0)
		return 1.0;

	double M[DET_SMALL_MAX_N * DET_SMALL_MAX_N];
	memcpy(M, A, sizeof(double) * (size_t)n * (size_t)n);

	double det = 1.0;

	for (uint32_t k = 0; k < n; k++) {
		uint32_t piv = k;
		double max_abs = fabs(M[k * n + k]);
		for (uint32_t i = k + 1; i < n; i++) {
			const double v = fabs(M[i * n + k]);
			if (v > max_abs) {
				max_abs = v;
				piv = i;
			}
		}
		if (max_abs == 0.0)
			return 0.0;

		if (piv != k) {
			for (uint32_t j = k; j < n; j++) {
				const double tmp = M[k * n + j];
				M[k * n + j] = M[piv * n + j];
				M[piv * n + j] = tmp;
			}
			det = -det;
		}

		const double pivot = M[k * n + k];
		det *= pivot;
		const double inv_pivot = 1.0 / pivot;

		for (uint32_t i = k + 1; i < n; i++) {
			const double f = M[i * n + k] * inv_pivot;
			if (f == 0.0)
				continue;
			for (uint32_t j = k + 1; j < n; j++)
				M[i * n + j] -= f * M[k * n + j];
		}
	}

	return det;
}
