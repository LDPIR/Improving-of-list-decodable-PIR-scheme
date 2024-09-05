// global.hpp

#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include "pch.hpp"
#include <string>  // 确保包含 <string> 头文件以使用 std::string

// 全局参数声明
extern int n;
extern int _index;
extern int m;
extern int w;
extern int t;
extern int k;
extern int b;
extern int h;
extern int D;

extern std::string E_name;
extern std::string data_name;
extern std::string label;

extern fmpz_t field_size_mpz;
extern fmpz_mod_ctx_t ctx;
extern fmpz_mod_mpoly_ctx_t mpoly_ctx;
extern fmpz_mod_mpoly_ctx_t database_ctx;

extern flint_rand_t state;
extern int** E;  // index set 
extern fmpz_t* DB;
extern const char** vars;
extern const char** database_vars;

#endif // GLOBAL_HPP
