// global.cpp

#include "global.hpp"

// 全局参数定义
int n = 1 << 16;
int _index = 2;
int m = 37;   
// 16 37
// 18 52
// 20 73
// 22 102
// 24 144
// 26 202


int w = 4;    
int t = 1;
int k = 12;
int b = 6;
int h = k - b;
int D = 2 * h - 1;  

std::string E_name =  "experiment_E16.txt";
std::string data_name = "experiment_database16.txt";
std::string label = "";

fmpz_t field_size_mpz;
fmpz_mod_ctx_t ctx;
fmpz_mod_mpoly_ctx_t mpoly_ctx;
fmpz_mod_mpoly_ctx_t database_ctx;

flint_rand_t state;
int** E;  // index set 
fmpz_t* DB;
const char** vars = new const char* [2];
const char** database_vars = new const char* [m];
