#include "2007_improving.hpp"
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>

fmpz_t p;

int _index = 2;
int n = 1 << 16;
int r = 1 << 8;
int s = 1 << 8;
std::string data_name = "experiment_database16.txt";
std::string label = "2";

fmpz_t field_size_mpz;
fmpz_mod_ctx_t ctx;
fmpz_mod_mpoly_ctx_t mpoly_ctx;

int t = 1;
int k = 8;
int b = 4;
int h = k - b;
int D = h - 1;  
// h > D
// b < k - floor{sqrt{kt}}
bool debug = true;

flint_rand_t state;
fmpz_t * DB;
const char ** vars = new const char *[2];


void Query(int n, int index, fmpz_t ** & Q, fmpz_t * & aux1, fmpz_t ** & aux2){
    // Q[0..k-1][0..r-1]
    // aux1[0..k-1]
    // aux2[0..t-1][0..r-1]
    int beta = int(index / s);
    // query the beta'th block instead of the index'th word  
    auto V = new fmpz_t * [t];
    auto lambda = new fmpz_t [k];


    for (int j = 0; j < k; j ++){
        fmpz_init(lambda[j]);
        // fmpz_randm(lambda[j], state, p);
        fmpz_set_ui(lambda[j], j + 1);
        fmpz_set(aux1[j], lambda[j]);
    }

    for (int s = 0; s < t; s ++){
        V[s] = new fmpz_t[r];
        for (int theta = 0; theta < r; theta ++){
            fmpz_init(V[s][theta]);
            fmpz_randm(V[s][theta] ,state, field_size_mpz);
            fmpz_set(aux2[s][theta], V[s][theta]);
        }
    }

    // Q_j = G(\lambda_j) = E(i) + sum_s lambda^s V_s
    fmpz_t temp;
    fmpz_init(temp);
    for (int j = 0; j < k; j ++){
        for (int theta = 0; theta < r; theta ++){
            if (theta == beta){
                fmpz_one(Q[j][theta]);
            }
            else{
                fmpz_zero(Q[j][theta]);
            }
            for (int s = t ; s > 0; s --){
                fmpz_pow_ui(temp, lambda[j], s);
                fmpz_mul(temp, temp, V[s - 1][theta]);
                fmpz_add(Q[j][theta], Q[j][theta], temp);
            }           
            fmpz_mod(Q[j][theta], Q[j][theta], field_size_mpz);  
        }
    }
    fmpz_clear(temp);
}

void Answer(int n, fmpz_t * DB, fmpz_t * _Q, fmpz_t * & _A){
    // _A[0..s-1]
    fmpz_t temp;
    fmpz_init(temp);
    for (int j = 0; j < s; j ++){
        fmpz_set_si(_A[j], 0);
    }
    for (int i = 0; i < r; i ++){
        for (int j = 0; j < s; j ++){
            fmpz_addmul(_A[j], _Q[i], DB[i * s + j]);
            fmpz_mod(_A[j],_A[j], field_size_mpz);            
        }
    }
}

void Reconstruct(int _index, fmpz_t ** _A, fmpz_t * aux1, fmpz_t ** aux2, std::vector<std::vector<FmpzWrapper>> & res){
    int count = 0;
    std::vector<H_struct> H;
    std::vector<H_struct> H_;

    H_struct temp_H;    
    for (int i = 0; i < k; i ++){
        temp_H.G.push_back(i);
    }
    H.push_back(temp_H);

    fmpz_t * R = new fmpz_t[k];
    for (int i = 0; i < k; i ++){
        fmpz_init(R[i]);
    }
    int HardRecover_count = 0;
    for (int c = 0; c < s; c ++){
        for (int i = 0; i < k; i ++){
            fmpz_mod(R[i], _A[i][c], field_size_mpz);
        }

        EasyRecover(H, H_, R, aux1); 

        // printf("round %d in ", c);

        if (H_.size() == 0){
            // printf("HardRecover\n");
            HardRecover(H, H_, R, aux1);
            HardRecover_count ++;
        }
        else{
            // printf("EasyRecover\n");
        }             
        H = H_;
        H_.clear();        
    }    
    
    printf("%d HardRecover, %d EasyRecover, %d rounds in total.\n", HardRecover_count, s - HardRecover_count, s);
    printf("%ld possible lists\n", H.size());
    
    for (int i = 0; i < H.size(); i ++){
        res.push_back(H[i].sigma);
    }    
}

// H_ is the output, should be empty before execution.
void EasyRecover(std::vector<H_struct> H, std::vector<H_struct> & H_, fmpz_t * R, fmpz_t * alpha){
    H_.clear();
    int _size = H.size();
    
    fmpz_mod_poly_t phi;
    fmpz_mod_poly_init(phi, ctx);
    fmpz_mod_poly_t temp_poly;
    fmpz_mod_poly_init(temp_poly, ctx);
    
    fmpz * temp_index; 
    temp_index = FLINT_ARRAY_ALLOC(t, fmpz);
    fmpz_t temp;
    fmpz_init(temp);

    for ( int i = 0; i < _size; i ++ ){
        // 用H[i].G的前t + 1个点所对应的index的做Lagrange Interpolation, phi(alpha_i)=R_i
        fmpz_mod_poly_zero(phi, ctx);
        // printf("interpolation points:\n");
        // for (int curr = 0; curr <= t; curr ++){
        //     printf("phi[");
        //     fmpz_print(alpha[H[i].G[curr]]);
        //     printf("] = ");
        //     fmpz_print(R[H[i].G[curr]]);
        //     printf("\n");
        // }

        for (int curr = 0; curr <= t; curr ++){
            // interpolate curr-th index
            for (int j = 0; j <= t; j ++){
                if (j < curr){
                    temp_index[j] = * alpha[H[i].G[j]];
                }
                if (j > curr){
                    temp_index[j - 1] = * alpha[H[i].G[j]];
                }
            }
            fmpz_mod_poly_product_roots_fmpz_vec(temp_poly, temp_index, t, ctx);
            
            fmpz_mod_poly_evaluate_fmpz(temp, temp_poly, alpha[H[i].G[curr]], ctx);
            
            fmpz_mod_poly_scalar_div_fmpz(temp_poly, temp_poly, temp, ctx);
            fmpz_mod_poly_scalar_addmul_fmpz(phi, temp_poly, R[H[i].G[curr]], ctx);
        }
        // assert the number of index such that phi(alpha_i) = R_i.
        int count = 0;
        
        std::vector<int> new_G;
        for (int curr = 0; curr < H[i].G.size(); curr ++){
            fmpz_mod_poly_evaluate_fmpz(temp, phi, alpha[H[i].G[curr]], ctx);
            if (fmpz_equal(temp, R[H[i].G[curr]])){
                new_G.push_back(H[i].G[curr]);
                count ++;
            }
        }    
        if ((count >= h) and (H[i].G.size() - count < h - t)){
            
            H_struct temp_H;
            temp_H.G = new_G;  
                
            temp_H.sigma = H[i].sigma;     
            fmpz_t temp;
            fmpz_init(temp);
            fmpz_t zero;
            fmpz_init(zero);
            fmpz_zero(zero);
            fmpz_mod_poly_evaluate_fmpz(temp, phi, zero, ctx);
            FmpzWrapper wrapper{temp};
            temp_H.sigma.push_back(wrapper);
                
            H_.push_back(temp_H);
        }
        else{
            H_.clear();
            return;
        }
    }
}

void HardRecover(std::vector<H_struct> H, std::vector<H_struct> & H_, fmpz_t * R, fmpz_t *aux1){
    H_.clear();

    fmpz_mod_poly_t * polys; 
    int count;
    Li_decoding(k, h, t, D, aux1, R, polys, count);
    
    // printf("output of Li-decoding:\n");
    // for (int i = 0; i < count; i ++){
    //     printf("poly[%d] = ", i);
    //     fmpz_mod_poly_print_pretty(polys[i], "X", ctx);
    //     printf("\n");
    // }
    
    
    
    for (int i = 0; i < count; i ++){
        // assert if poly[i] satisfy at least h points (in each G of H).
        int correct_point_num = 0;
        for (int j = 0; j < H.size(); j ++){
            std::vector<int> new_G;
            for (int _k = 0; _k < H[j].G.size(); _k ++){
                fmpz_t temp;
                fmpz_init(temp);
                fmpz_mod_poly_evaluate_fmpz(temp, polys[i], aux1[H[j].G[_k]], ctx);
                if (fmpz_cmp(temp, R[H[j].G[_k]]) == 0){
                    new_G.push_back(H[j].G[_k]);
                }
            }
            if (new_G.size() >= h){
                H_struct temp_H;
                temp_H.G = new_G;  
                
                temp_H.sigma = H[j].sigma;     
                fmpz_t temp;
                fmpz_init(temp);
                fmpz_t zero;
                fmpz_init(zero);
                fmpz_zero(zero);
                fmpz_mod_poly_evaluate_fmpz(temp, polys[i], zero, ctx);
                FmpzWrapper wrapper{temp};
                temp_H.sigma.push_back(wrapper);
                
                H_.push_back(temp_H);
            }
        }
    }    
}

void Li_decoding(int total_num, int corr_num, int degree, int D, fmpz_t * lambda, fmpz_t * alpha, fmpz_mod_poly_t * & polys, int & count){    

    int l = int(D / degree);

// 
    for (int i = 0; i < total_num; i ++){
        fmpz_mod(lambda[i], lambda[i], field_size_mpz);
        fmpz_mod(alpha[i], alpha[i], field_size_mpz);
        // fmpz_mod(beta[i], beta[i], field_size_mpz);
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
    fmpz_mod_mat_init(A, total_num + 1, num_coeff, ctx); 
    fmpz_mod_mat_init(X, num_coeff, 1, ctx);
    fmpz_mod_mat_init(B, total_num + 1, 1, ctx);
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
                fmpz_zero(valx);
                fmpz_zero(valy);
                fmpz_mod_pow_ui(valx, lambda[_r], deg_x, ctx);
                fmpz_mod_pow_ui(valy, alpha[_r], deg_y, ctx);
                // printf("deg_x = %d, deg_y = %d, ", deg_x, deg_y);
                // printf("val_x = ");                
                // fmpz_print(valx);
                // printf(" val_y = ");                
                // fmpz_print(valy);
                // printf("\n");
                // fmpz_mod_mul(val, valx, valy, ctx);                           

                // printf(" val = ");
                // fmpz_print(val);
                // printf("\n");
                fmpz_mod_mat_set_entry(A, _r, _c, val, ctx);
                _c = _c + 1;
            }    
        }
    }
    // add one row under original matrix A to avoid output zero vector x.
    fmpz_t one;
    fmpz_init(one);
    fmpz_set_ui(one, 1);
    for (int _c = 0; _c < num_coeff; _c ++){
        fmpz_mod_mat_set_entry(A, total_num, _c, one, ctx);
    }
    fmpz_mod_mat_set_entry(B, total_num, 0, one, ctx);

    int res = fmpz_mod_mat_can_solve(X, A, B, ctx);
    
    printf("A = ");
    fmpz_mod_mat_print_pretty(A, ctx);
    printf("\n");
    printf("X = ");
    fmpz_mod_mat_print_pretty(X, ctx);
    printf("\n");
    printf("B = ");
    fmpz_mod_mat_print_pretty(B, ctx);
    printf("\n");


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
    // printf("\n Q_poly:\n");
    // fmpz_mod_mpoly_print_pretty(Q_poly, vars, mpoly_ctx);
    // printf("\n");

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
                // fmpz_t neg;
                // fmpz_init(neg);
                // fmpz_set(neg, field_size_mpz);
                // fmpz_mod_mul(div, div, neg, ctx);
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


int main(){
    fmpz_init(field_size_mpz);
    fmpz_set_d_2exp(field_size_mpz, (double)1, 128);
    fmpz_nextprime(field_size_mpz, field_size_mpz, 1);
    fmpz_mod_ctx_init(ctx, field_size_mpz);
    fmpz_mod_mpoly_ctx_init(mpoly_ctx, 2, ORD_LEX, field_size_mpz);

    DB = new fmpz_t[r * s];
    flint_randinit(state);

    std::string data_path = DATA_DIR;
    data_path += data_name;  

    std::ifstream DBifs(data_path);
    // 检查文件是否成功打开
    if (!DBifs.is_open()) {
        std::cerr << "Error: Could not open file " << data_name << std::endl;
        return 1; // 返回错误代码
    }

    // 读取数据库DB
    int hi;
    int lo;

    for (int i = 0; i < n; i ++){
        DBifs >> hi;
        DBifs >> lo;
        fmpz_set_uiui(DB[i], hi, lo);
    }

    DBifs.close(); 

    for (int i = n; i < r * s; i ++){
        fmpz_zero(DB[i]);
    }

    vars[0] = "X";
    vars[1] = "Z";


    fmpz_t * aux1 = new fmpz_t[k];
    for (int j = 0; j < k; j ++){
        fmpz_init(aux1[j]);
    }
    fmpz_t ** aux2 = new fmpz_t *[t];
    for (int s = 0; s < t; s ++){
        aux2[s] = new fmpz_t[r];
        for (int j = 0; j < r; j ++){
            fmpz_init(aux2[s][j]);
        }
    }


    fmpz_t ** Q = new fmpz_t *[k];
    for (int j = 0; j < k; j ++){
        Q[j] = new fmpz_t[r];
    }

    printf("before query\n");

    auto qu_start = clock();
    Query(n, _index, Q, aux1, aux2);
    auto qu_end = clock();

    printf("after query\n");

    fmpz_t ** _A = new fmpz_t *[k];
    for (int j = 0 ; j < k; j ++){
        _A[j] = new fmpz_t [s];
        for (int i = 0; i < s; i ++){
            fmpz_init(_A[j][i]);   
        }
    }

    printf("before answer\n");
    auto ans_start = clock();    
    for (int j = 0; j < k; j ++){
        Answer(n, DB, Q[j], _A[j]);
    }
    auto ans_end = clock();
    printf("after answer\n");

    // generate t random numbers in [0..k-1] to simulate incorrect answers
    std::vector<int> random_numbers = generate_unique_random_numbers(b, k);
    fmpz_t addition;
    fmpz_init(addition);
    for (int num : random_numbers) {
    // for (int num = 0; num < b; num ++){
        for (int j = 0; j < s; j ++){
            // fmpz_set_ui(addition, 2);
            // fmpz_add(_A[num][j], _A[num][j], addition);
            fmpz_randm(_A[num][j], state, field_size_mpz);
        }
    }
    
    std::vector<std::vector<FmpzWrapper>> res;

    printf("before reconstruction\n");
    auto re_start = clock();
    Reconstruct(_index, _A , aux1, aux2, res);
    auto re_end = clock();
    printf("after reconstruction\n");

    printf("(t,k,b) = (%d,%d,%d)\n", t, k, b);
    printf("database size = 2^%d blocks\n", __builtin_ctz(n));
    printf("blcok size = 128 bits\n");
    printf("query takes time %f s\n", (double)(qu_end - qu_start)/CLOCKS_PER_SEC);
    printf("answer takes time %f s\n", (double)(ans_end - ans_start)/CLOCKS_PER_SEC);
    printf("reconstruct takes time %f s\n", (double)(re_end - re_start)/CLOCKS_PER_SEC);


    // fmpz_t ** Q[0..k-1][0..r-1];
    // fmpz_t ** _A[0..k-1][0..s-1];
    
    std::string output_query = "output_query" + label + data_name; 
    std::ofstream outFile(output_query);
    if (outFile.is_open()) {
        for (int j = 0; j < k; j ++){
            for (int theta = 0; theta < r; theta ++){
                char* str = fmpz_get_str(NULL, 10, Q[j][theta]); // 将 fmpz_t 转换为字符串
                outFile << str << std::endl;
                flint_free(str); // 释放由 fmpz_get_str 分配的内存
            }
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file";
    }

    std::string output_answer = "output_answer" + label + data_name; 
    std::ofstream outFile2(output_answer);
    if (outFile2.is_open()) {
        for (int j = 0; j < k; j ++){
            for (int i = 0; i < s; i ++){
                char* str = fmpz_get_str(NULL, 10, _A[j][i]); // 将 fmpz_t 转换为字符串
                outFile2 << str << std::endl;
                flint_free(str); // 释放由 fmpz_get_str 分配的内存
            }
        }
        outFile2.close();
    } else {
        std::cerr << "Unable to open file";
    }

    flint_randclear(state);
    return 0;
}


