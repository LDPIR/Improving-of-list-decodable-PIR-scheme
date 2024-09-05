#include "pch.hpp"
#include "global.hpp"

void ServerPreprocess(fmpz_t * DB, int ** E, fmpz_mod_mpoly_t & F, fmpz_mod_mpoly_t * & partialF);

void Answer_With_Preprocess(int n, fmpz_mod_mpoly_t F, fmpz_mod_mpoly_t * partialF, fmpz_t * _Q, fmpz_t & A_val, fmpz_t *& A_partial);

void Answer_Without_Preprocess(int n, fmpz_t * DB, fmpz_t * _Q, fmpz_t & A_val, fmpz_t *& A_partial);
