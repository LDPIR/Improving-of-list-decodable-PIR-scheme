#include "pch.hpp"
#include "global.hpp"


void Reconstruct(int _index, fmpz_t * A_vals, fmpz_t ** A_partials, fmpz_t * aux1, fmpz_t ** aux2, fmpz_t * & res);

void Li_decoding(int total_num, int corr_num, int degree, int D, fmpz_t * lambda, fmpz_t * alpha, fmpz_t * beta, fmpz_mod_poly_t * & polys, int & count);
