#include "query.hpp"

// aux1 = lambda_1,...,lambda_k   aux2 = V_1,..,V_t
void Query(int n, int _index, fmpz_t ** & Q, fmpz_t * & aux1, fmpz_t ** & aux2){
    // Q[0..k-1][0..m-1]
    // aux1[0..k-1]
    // aux2[0..t-1][0..m-1]
    auto V = new fmpz_t * [t];
    auto lambda = new fmpz_t [k];
    for (int s = 0; s < t; s ++){
        V[s] = new fmpz_t[m];
        for (int theta = 0; theta < m; theta ++){
            fmpz_init(V[s][theta]);
            fmpz_randm(V[s][theta] ,state, field_size_mpz);
            fmpz_set(aux2[s][theta], V[s][theta]);
        }
    }
    for (int j = 0; j < k; j ++){
        fmpz_init(lambda[j]);
        fmpz_randm(lambda[j], state, field_size_mpz);
        fmpz_init(aux1[j]);
        fmpz_set(aux1[j], lambda[j]);
    }
    
    // Q_j = G(\lambda_j) = E(i) + sum_s lambda^s V_s
    fmpz_mod_poly_t *G =new fmpz_mod_poly_t [m];
    for (int theta = 0; theta < m; theta ++){
        fmpz_mod_poly_init(G[theta], ctx); 
    }

    for (int deg = 0; deg <= t; deg ++){
        if (deg == 0){
            for (int theta = 0; theta < m; theta ++){
                fmpz_mod_poly_zero(G[theta], ctx);
            }
            for (int _w = 0; _w < w; _w ++){
                fmpz_mod_poly_set_coeff_ui(G[E[_index][_w]], deg, 1, ctx);            
            }
        }   
        else {
            for (int theta = 0; theta < m; theta ++){
                fmpz_mod_poly_set_coeff_fmpz(G[theta], deg, V[deg - 1][theta], ctx);
            }
        }
    }

    // for (int theta = 0; theta < m; theta ++){
    //     fmpz_mod_poly_print_pretty(G[theta], "X", ctx); 
    //     printf("\n");
    // }

    for (int j = 0 ; j < k; j ++){
        for (int theta = 0; theta < m; theta ++){
            fmpz_mod_poly_evaluate_fmpz(Q[j][theta], G[theta], lambda[j], ctx);
        }
    }
}
