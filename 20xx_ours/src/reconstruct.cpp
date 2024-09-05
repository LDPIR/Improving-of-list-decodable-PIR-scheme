#include "reconstruct.hpp"

void Reconstruct(int _index, fmpz_t * A_vals, fmpz_t ** A_partials, fmpz_t * aux1, fmpz_t ** aux2, fmpz_t * & res){
    fmpz_t temp;
    fmpz_t * temp_vec = new fmpz_t[m];
    fmpz_init(temp);
    for (int theta = 0; theta < m; theta ++){
        fmpz_init(temp_vec[theta]);
    }

    fmpz_t * alpha = new fmpz_t[k];
    fmpz_t * beta = new fmpz_t[k];
    for (int j = 0; j < k; j ++){
        fmpz_init(alpha[j]);
        fmpz_init(beta[j]);
    }

    for (int j = 0; j < k; j ++){
        // printf(" j = %d, \n", j);
        fmpz_set(alpha[j], A_vals[j]);
        fmpz_set_si(beta[j], 0);
        fmpz_set_si(temp, 0);
        
        for (int s = 0; s < t; s ++){
            fmpz_pow_ui(temp, aux1[j], s);
            fmpz_mul_ui(temp, temp, s + 1);
            for (int theta = 0; theta < m; theta ++){
                fmpz_mul(temp_vec[theta], temp, aux2[s][theta]);
            }
        }
        for (int theta = 0; theta < m; theta ++){
            fmpz_addmul(beta[j], temp_vec[theta], A_partials[j][theta]);
        }
    }    

    // De_decoding();
    fmpz_mod_poly_t * polys;
    int count;

    printf("before list decoding\n");
    Li_decoding(k, k - b, w, D, aux1, alpha, beta, polys, count);
    printf("after list decoding\n");



    for (int s = 0; s < count; s ++){
        printf("polys[%d] = ", s);
        fmpz_mod_poly_print_pretty(polys[s], "X", ctx);
        printf("\n");
    }

}

void Li_decoding(int total_num, int corr_num, int degree, int D, fmpz_t * lambda, fmpz_t * alpha, fmpz_t * beta, fmpz_mod_poly_t * & polys, int & count){    

    int l = int(D / degree);

    for (int i = 0; i < total_num; i ++){
        fmpz_mod(lambda[i], lambda[i], field_size_mpz);
        fmpz_mod(alpha[i], alpha[i], field_size_mpz);
        fmpz_mod(beta[i], beta[i], field_size_mpz);
    }
// 

    // printf("D = %d, degree = %d, l = %d \n", D, degree, l);

    // printf("lambda = ");
    // for (int i = 0; i < total_num; i ++){
    //     fmpz_print(lambda[i]);
    //     printf(" ");
    // }
    // printf("\n");

    // printf("alpha = ");
    // for (int i = 0; i < total_num; i ++){
    //     fmpz_print(alpha[i]);
    //     printf(" ");
    // }
    // printf("\n");

    // printf("beta = ");
    // for (int i = 0; i < total_num; i ++){
    //     fmpz_print(beta[i]);
    //     printf(" ");
    // }
    // printf("\n");


    int num_coeff = 0;
    for (int deg_y = 0; deg_y <= l; deg_y ++){
        for (int deg_x = 0; deg_x <= D - deg_y * degree; deg_x ++ ){
            num_coeff = num_coeff + 1;
        }    
    }
    fmpz_t * Q = new fmpz_t [num_coeff];
    for (int j = 0; j < num_coeff; j ++){
        fmpz_init(Q[j]);
    }
    
    // solve Ax = 0 where x is the coeffs of Q
    // To avoid output a zero vector x, add one row under origin matrix A.

    fmpz_mod_mat_t A;
    fmpz_mod_mat_t X;
    fmpz_mod_mat_t B;
    fmpz_mod_mat_init(A, 2 * total_num + 1, num_coeff, ctx); 
    fmpz_mod_mat_init(X, num_coeff, 1, ctx);
    fmpz_mod_mat_init(B, 2 * total_num + 1, 1, ctx);
    fmpz_mod_mat_zero(B, ctx);

    // set the first n rows with vals
    for (int _r = 0; _r < total_num; _r ++){
        int _c = 0;
        fmpz_t valx;
        fmpz_t valy;
        fmpz_t val;
        fmpz_init(val);
        fmpz_init(valx);
        fmpz_init(valy);
        
        for (int deg_y = 0; deg_y <= l; deg_y ++){
            for (int deg_x = 0; deg_x <= D - deg_y * degree; deg_x ++){
                fmpz_mod_pow_ui(valx, lambda[_r], deg_x, ctx);
                fmpz_mod_pow_ui(valy, alpha[_r], deg_y, ctx);
                // printf("deg_x = %d, deg_y = %d, ", deg_x, deg_y);
                // printf("val_x = ");                
                // fmpz_print(valx);
                // printf(" val_y = ");                
                // fmpz_print(valy);

                fmpz_mod_mul(val, valx, valy, ctx);                           

                // printf(" val = ");
                // fmpz_print(val);
                // printf("\n");

                fmpz_mod_mat_set_entry(A, _r, _c, val, ctx);
                _c = _c + 1;
            }    
        }
    }
    // set the last n rows with diffs
    for (int _r = 0; _r < total_num; _r ++){
        int _c = 0;
        fmpz_t valx;
        fmpz_t valy;
        fmpz_t val;
        fmpz_init(val);
        fmpz_init(valx);
        fmpz_init(valy);
        
        for (int deg_y = 0; deg_y <= l; deg_y ++){
            for (int deg_x = 0; deg_x <= D - deg_y * degree; deg_x ++){
                // diff of x
                if (deg_x > 0){
                    fmpz_mod_pow_ui(valx, lambda[_r], deg_x - 1, ctx);
                    fmpz_mod_mul_ui(valx, valx, deg_x, ctx);
                }
                else{
                    fmpz_set_ui(valx, 0);
                }
                fmpz_mod_pow_ui(valy, alpha[_r], deg_y, ctx);
                fmpz_mod_mul(val, valx, valy, ctx);                           
                // diff of y
                fmpz_mod_pow_ui(valx, lambda[_r], deg_x , ctx);
                if (deg_y > 0){
                    fmpz_mod_pow_ui(valy, alpha[_r], deg_y - 1, ctx);
                    fmpz_mod_mul_ui(valy, valy, deg_y, ctx);
                    fmpz_mod_mul(valy, valy, beta[_r], ctx);
                }
                fmpz_mod_mul(valy, valx, valy, ctx);
                
                fmpz_mod_add(val, val, valy, ctx);
                fmpz_mod_mat_set_entry(A, _r + total_num, _c, val, ctx);
                _c = _c + 1;
            }    
        }
    }

    // add one row under original matrix A to avoid output zero vector x.
    fmpz_t one;
    fmpz_init(one);
    fmpz_set_ui(one, 1);
    for (int _c = 0; _c < num_coeff; _c ++){
        fmpz_mod_mat_set_entry(A, 2 * total_num, _c, one, ctx);
    }
    fmpz_mod_mat_set_entry(B, 2 * total_num, 0, one, ctx);

    int res = fmpz_mod_mat_can_solve(X, A, B, ctx);
    
    // printf("A = ");
    // fmpz_mod_mat_print_pretty(A, ctx);
    // printf("\n");
    // printf("X = ");
    // fmpz_mod_mat_print_pretty(X, ctx);
    // printf("\n");
    // printf("B = ");
    // fmpz_mod_mat_print_pretty(B, ctx);
    // printf("\n");



    if (res == 0){
        printf("\n\n\n cannot solve \n\n\n");
        return;
    }
        
    // fmpz_mod_mat_print_pretty(A, ctx);
    // fmpz_mod_mat_print_pretty(B, ctx);
    // fmpz_mod_mat_print_pretty(X, ctx); 
    fmpz_mod_mpoly_t Q_poly;    
    fmpz_mod_mpoly_init(Q_poly, mpoly_ctx);
    int _c = 0;
    ulong * exp = new ulong [2];    
    for (int deg_y = 0; deg_y <= l; deg_y ++){
        for (int deg_x = 0; deg_x <= D - deg_y * degree; deg_x ++){ 
            exp[0] = deg_x;
            exp[1] = deg_y;
            fmpz_mod_mpoly_set_coeff_fmpz_ui(Q_poly, fmpz_mod_mat_entry(X, _c, 0), exp, mpoly_ctx);
            _c = _c + 1;
        }    
    }
    printf("\n Q_poly:\n");
    fmpz_mod_mpoly_print_pretty(Q_poly, vars, mpoly_ctx);
    printf("\n");

    fmpz_mod_mpoly_factor_t factor;
    fmpz_mod_mpoly_factor_init(factor, mpoly_ctx);
    fmpz_mod_mpoly_factor(factor, Q_poly, mpoly_ctx);

    fmpz_t temp;
    fmpz_mod_mpoly_t temp_mpoly;
    fmpz_init(temp);
    fmpz_mod_mpoly_init(temp_mpoly, mpoly_ctx);
    slong length = fmpz_mod_mpoly_factor_length(factor, mpoly_ctx);
    fmpz_mod_mpoly_factor_get_constant_fmpz(temp, factor, mpoly_ctx);
    count = length;
    // printf("factor of Q_poly:\n");
    // fmpz_print(temp);
    // printf("\n");
    // for(slong i = 0; i < length; i ++){
    //     fmpz_mod_mpoly_factor_get_base(temp_mpoly, factor, i, mpoly_ctx);        
    //     fmpz_mod_mpoly_print_pretty(temp_mpoly, vars, mpoly_ctx);
    //     printf("\n");

    //     // printf("degree of X: %ld\n", fmpz_mod_mpoly_degree_si(temp_mpoly, 0, mpoly_ctx));
    //     // printf("degree of Z: %ld\n", fmpz_mod_mpoly_degree_si(temp_mpoly, 1, mpoly_ctx));
    // }

    polys = new fmpz_mod_poly_t [length];
    for (int i = 0; i < length; i ++){
        fmpz_mod_poly_init(polys[i], ctx);
        fmpz_mod_mpoly_factor_get_base(temp_mpoly, factor, i,  mpoly_ctx);

        // assert if the factor is in form of Z+f(X)
        if (fmpz_mod_mpoly_degree_si(temp_mpoly, 1, mpoly_ctx) == 1 ){
            int deg_x = fmpz_mod_mpoly_degree_si(temp_mpoly, 0, mpoly_ctx);
            bool flag = true;
            fmpz_t temp;
            fmpz_init(temp);

            ulong deg[2];
            deg[1] = 1;
            for (ulong _deg = 1; _deg < deg_x; _deg ++){
                deg[0] = _deg;
                fmpz_mod_mpoly_get_coeff_fmpz_ui(temp, temp_mpoly, deg, mpoly_ctx);
                if (not fmpz_is_zero(temp)){
                    flag = false;
                }
            }            

            if (flag){
                fmpz_t temp;
                fmpz_init(temp);
                fmpz_t div;
                fmpz_t temp2;
                fmpz_init(div);
                fmpz_init(temp2);

                deg[0] = 0;
                deg[1] = 1;
                fmpz_mod_mpoly_get_coeff_fmpz_ui(div, temp_mpoly, deg, mpoly_ctx);
                
                deg[1] = 0;
                fmpz_sub(div, field_size_mpz, div);


                for (ulong _deg = 0; _deg <= deg_x; _deg ++){
                    deg[0] = _deg;
                    fmpz_mod_mpoly_get_coeff_fmpz_ui(temp, temp_mpoly, deg, mpoly_ctx); 
                    fmpz_mod_divides(temp2, temp, div, ctx);

                    fmpz_mod_poly_set_coeff_fmpz(polys[i], _deg, temp2, ctx);                                   
                }
            }
        }
    }
}

