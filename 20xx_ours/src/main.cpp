#include "main.hpp"
#include "global.hpp"

int main(){
    fmpz_init(field_size_mpz);
    fmpz_set_d_2exp(field_size_mpz, (double)1, 128);
    fmpz_nextprime(field_size_mpz, field_size_mpz, 1);
    // fmpz_set_ui(field_size_mpz, 97);
    fmpz_mod_ctx_init(ctx, field_size_mpz);
    fmpz_mod_mpoly_ctx_init(mpoly_ctx, 2, ORD_LEX, field_size_mpz);
    fmpz_mod_mpoly_ctx_init(database_ctx, m, ORD_LEX, field_size_mpz);

    E = new int *[n];
    DB = new fmpz_t[n];
    flint_randinit(state);
    for (int i = 0; i < n; i ++){
        E[i] = new int[w];
        fmpz_init(DB[i]);
    }
    
    std::string E_path = DATA_DIR;
    E_path += E_name;  
    std::string data_path = DATA_DIR;
    data_path += data_name;  


    // 打开 "data_E" 文件
    std::ifstream Eifs(E_path);
    // 检查文件是否成功打开
    if (!Eifs.is_open()) {
        std::cerr << "Error: Could not open file" << E_name << std::endl;
        // printf("E_path = %s\n", E_path);
        return 1; // 返回错误代码
    }
    // 读取文件
    // 读取Mapping E
    for (int i = 0; i < n; i ++){
        for (int _w = 0; _w < w; _w ++){
            Eifs >> E[i][_w];
        }
    }
    Eifs.close(); 

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



    vars[0] = "X";
    vars[1] = "Z";

    for (int theta = 0; theta < m; theta ++){
        char buffer[20];
        std::sprintf(buffer, "z_%d ", theta);
        database_vars[theta] = new char[strlen(buffer) + 1];
        std::strcpy(const_cast<char *>(database_vars[theta]), buffer);
    }


    fmpz_t * aux1 = new fmpz_t[k];
    fmpz_t ** aux2 = new fmpz_t *[t];
    for (int s = 0; s < t; s ++){
        aux2[s] = new fmpz_t[m];
    }

    fmpz_t ** Q = new fmpz_t *[k];
    for (int j = 0; j < k; j ++){
        Q[j] = new fmpz_t[m];
        for (int theta = 0; theta < m; theta ++){
            fmpz_init(Q[j][theta]);
        }
    }

    printf("start query\n");

    auto qu_start = clock();    
    Query(n, _index, Q, aux1, aux2);
    auto qu_end = clock();    
    printf("end query\n");





    printf("start ans\n");

    fmpz_t * A_val = new fmpz_t[k];
    fmpz_t ** A_partial = new fmpz_t *[k];
    for (int j = 0 ; j < k; j ++){
        fmpz_init(A_val[j]);
        A_partial[j] = new fmpz_t[m];
        for (int theta = 0; theta < m; theta ++){
            fmpz_init(A_partial[j][theta]);
        }
    }


    auto ans_without_pre_start = clock();    
    for (int j = 0; j < k; j ++){
        printf("j  = %d\n", j);
        Answer_Without_Preprocess(n, DB, Q[j], A_val[j], A_partial[j]);
    }
    auto ans_without_pre_end = clock();    

    fmpz_mod_mpoly_t F;
    fmpz_mod_mpoly_t * partialF = new fmpz_mod_mpoly_t [m];





    // auto pre_start = clock();    
    // ServerPreprocess(DB, E, F, partialF);
    // auto pre_end = clock();    


    // auto ans_with_pre_start = clock();    
    // for (int j = 0; j < k; j ++){
    //     Answer_With_Preprocess(n, F, partialF, Q[j], A_val_new[j], A_partial_new[j]);
    // }
    // auto ans_with_pre_end = clock();    

    printf("end ans\n");


    // bool flag = true;
    // for (int j = 0 ; j < k; j ++){
    //     if (fmpz_cmp(A_val[j], A_val_new[j]) != 0){
    //         flag = false;
    //     }
        
    //     for (int theta = 0; theta < m; theta ++){
    //         if (fmpz_cmp(A_partial[j][theta], A_partial_new[j][theta]) != 0){
    //             flag = false;
    //         }
    //     }
    //     if  (flag == false){
    //         break;
    //     }
    // }
    // printf("equivity flag = %s\n", flag ? "true" : "false");



    printf("start reconstruct\n");
    auto re_start = clock();    


    fmpz_t * A_val_new = new fmpz_t[k];
    fmpz_t ** A_partial_new = new fmpz_t *[k];
    for (int j = 0 ; j < k; j ++){
        fmpz_init(A_val_new[j]);
        A_partial_new[j] = new fmpz_t[m];
        for (int theta = 0; theta < m; theta ++){
            fmpz_init(A_partial_new[j][theta]);
        }
    }



    fmpz_t * res = new fmpz_t[100];
    for (int i = 0; i < 100; i ++){
        fmpz_init(res[i]);
    }    

    auto re_mid = clock();    
    Reconstruct(_index, A_val, A_partial, aux1, aux2, res);
    auto re_end = clock();    

    printf("end reconstruct\n");



    printf("(t,k,b) = (%d,%d,%d)\n", t, k, b);
    printf("database size = 2^%d blocks\n", __builtin_ctz(n));
    printf("blcok size = 128 bits\n");
    printf("query takes time %f s\n", (double)(qu_end - qu_start)/CLOCKS_PER_SEC);
    printf("answer without pre takes time %f s\n", (double)(ans_without_pre_end - ans_without_pre_start)/CLOCKS_PER_SEC);
    
    printf("reconstruct1 takes time %f s\n", (double)(re_mid - re_start)/CLOCKS_PER_SEC);
    printf("reconstruct2 takes time %f s\n", (double)(re_end - re_mid)/CLOCKS_PER_SEC);




    // printf("i = %d \n", _index);
    // printf("E[i]:\n");
    // for (int i = _index; i < _index + 1; i ++){
    //     for (int _w = 0; _w < w; _w ++){
    //         printf("%d ", E[i][_w]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // printf("aux1(lambda)：\n");
    // for (int j = 0; j < k; j ++){
    //     fmpz_print(aux1[j]);
    //     printf(" ");
    // }
    // printf("\n");

    // printf("aux2(V)：\n");
    // for (int s = 0; s < t; s ++){
    //     for (int theta = 0; theta < m; theta ++){
    //         fmpz_print(aux2[s][theta]);
    //         printf(" ");
    //     }
    //     printf("\n");
    // }
    // printf("\n");


    // printf("Q：\n");
    // for (int j = 0; j < k; j ++){
    //     for (int theta = 0; theta < m; theta ++){
    //         fmpz_print(Q[j][theta]);
    //         printf(" ");
    //     }
    //     printf("\n");
    // }
    // printf("Answer:\n");
    // for (int j = 0; j < k; j ++){
    //     if (j == 1){
    //     printf("A_val[%d] = ",j);
    //     fmpz_print(A_val[j]);
    //     printf("\n");
        
    //     printf("A_partial[%d] = ",j);
    //     for (int theta = 0; theta < m; theta ++){
    //         fmpz_print(A_partial[j][theta]);
    //         printf(" ");
    //     }
    //     printf("\n");
    //     }
    // }
    
    std::string query_file = "output_query" + label + data_name + ".txt";
    std::ofstream outFile(query_file);
    if (outFile.is_open()) {
        for (int j = 0; j < k; j ++){
            for (int theta = 0; theta < m; theta ++){
                char* str = fmpz_get_str(NULL, 10, Q[j][theta]); // 将 fmpz_t 转换为字符串
                outFile << str << std::endl;
                flint_free(str); // 释放由 fmpz_get_str 分配的内存
            }
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file";
    }

    std::string answer_file = "output_answer" + label + data_name + ".txt";
    std::ofstream outFile2(answer_file);
    if (outFile2.is_open()) {
        for (int j = 0; j < k; j ++){
            char* str = fmpz_get_str(NULL, 10, A_val[j]); 
            outFile2 << str << std::endl;
            flint_free(str); 

            for (int theta = 0; theta < m; theta ++){
                char* str = fmpz_get_str(NULL, 10, A_partial[j][theta]); 
                outFile2 << str << std::endl;
                flint_free(str); 
            }
        }
        outFile2.close();
    } else {
        std::cerr << "Unable to open file";
    }


    flint_randclear(state);
    for (auto i = 1; i < n; i ++){
        delete[] E[i];
    }
    delete[] E;
    
    return 0;
}


