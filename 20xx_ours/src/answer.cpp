#include "answer.hpp"

void ServerPreprocess(fmpz_t * DB, int ** E, fmpz_mod_mpoly_t & F, fmpz_mod_mpoly_t * & partialF){
        
    fmpz_mod_mpoly_init(F, database_ctx);
    for (int theta = 0; theta < m; theta ++){
        fmpz_mod_mpoly_init(partialF[theta], database_ctx);
    }

    // construct polynomial F
    ulong exp[m];
    for (int i = 0; i < n; i ++){
        for (int theta = 0; theta < m; theta ++){
            exp[theta] = 0;
        }
        for (int _w = 0; _w < w; _w ++){
            int temp_index = E[i][_w];
            exp[temp_index] = 1;                    
        }
        fmpz_mod_mpoly_set_coeff_fmpz_ui(F, DB[i], exp, database_ctx);
    }

    // derivate 
    for (int theta = 0; theta < m; theta ++){
        fmpz_mod_mpoly_derivative(partialF[theta], F, theta, database_ctx);
    }
}


void Answer_With_Preprocess(int n, fmpz_mod_mpoly_t F, fmpz_mod_mpoly_t * partialF, fmpz_t * _Q, fmpz_t & A_val, fmpz_t *& A_partial) {
    // 分配 evaluate_point 数组
    fmpz ** evaluate_point = new fmpz *[m];
    
    // 初始化 evaluate_point 数组中的每个 fmpz 元素
    for (int theta = 0; theta < m; theta++) {
        evaluate_point[theta] = new fmpz;
        fmpz_init(evaluate_point[theta]);
        fmpz_set(evaluate_point[theta], _Q[theta]);
    }
    
    // 设置 evaluate_point_pointer
    fmpz * const * evaluate_point_pointer = evaluate_point;

    // 初始化 A_val
    fmpz_init(A_val);
    
    // 调用函数
    fmpz_mod_mpoly_evaluate_all_fmpz(A_val, F, evaluate_point_pointer, database_ctx);
    for (int theta = 0; theta < m; theta ++){
        fmpz_mod_mpoly_evaluate_all_fmpz(A_partial[theta], partialF[theta], evaluate_point_pointer, database_ctx);    
    }

    // 清理内存
    for (int theta = 0; theta < m; theta++) {
        fmpz_clear(evaluate_point[theta]);
        delete evaluate_point[theta];
    }
    delete[] evaluate_point;
}

void Answer_Without_Preprocess(int n, fmpz_t * DB, fmpz_t * _Q, fmpz_t & A_val, fmpz_t *& A_partial){

    // printf("query = ");
    // for (int theta = 0; theta < m; theta ++){
    //     fmpz_print(_Q[theta]);
    //     printf(" ");
    // }
    // printf("\n");

    fmpz_t temp;
    fmpz_init(temp);

    fmpz_set_si(A_val, 0);
    for (int theta = 0; theta < m; theta ++){
        fmpz_set_si(A_partial[theta], 0);
    }    
    
    for (int i = 0; i < n; i ++){
        fmpz_set(temp, DB[i]);
        // A_val
        for (int _w = 0; _w < w; _w ++){
            fmpz_mul(temp, temp, _Q[E[i][_w]]);
            fmpz_mod(temp, temp, field_size_mpz);
        }
        fmpz_add(A_val, A_val, temp);
        fmpz_mod(A_val, A_val, field_size_mpz);

        // A_diff
        for (int _w = 0; _w < w; _w ++){
            fmpz_set(temp, DB[i]);
            for (int mul_index = 0; mul_index < w; mul_index ++){
                if (mul_index != _w){
                    fmpz_mod_mul(temp, temp, _Q[E[i][mul_index]], ctx);
                }
            }
            fmpz_mod_add(A_partial[E[i][_w]], A_partial[E[i][_w]], temp, ctx);
        }
    }    
}
