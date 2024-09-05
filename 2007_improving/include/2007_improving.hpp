#include <flint/fmpz_mod_types.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/fmpz_mod_mpoly_factor.h>

#include <gmp.h>

#include <time.h>
#include <thread>
class FmpzWrapper{
    public:
    fmpz_t value;

    FmpzWrapper(){
        fmpz_init(value);
    }

    FmpzWrapper(const fmpz_t& other){
        fmpz_init(value);
        fmpz_set(value, other);
    }


    FmpzWrapper(const FmpzWrapper& other){
        fmpz_init(value);
        fmpz_set(value, other.value);
    }

    FmpzWrapper& operator=(const FmpzWrapper & other){
        if (this != & other){
            fmpz_set(value, other.value);
        }
        return * this;
    }

    ~FmpzWrapper(){
        fmpz_clear(value);
    }
};


struct H_struct {
    std::vector<int> G;
    std::vector<FmpzWrapper> sigma;
};


void EasyRecover(std::vector<H_struct> H, std::vector<H_struct> & H_, fmpz_t * R, fmpz_t * aux1);                

void HardRecover(std::vector<H_struct> H, std::vector<H_struct> & H_, fmpz_t * R, fmpz_t * aux1);                

void Li_decoding(int total_num, int corr_num, int degree, int D, fmpz_t * lambda, fmpz_t * alpha, fmpz_mod_poly_t * & polys, int & count);


std::vector<int> generate_unique_random_numbers(int t, int k) {
    // 检查 t 和 k 的合法性
    if (t > k) {
        throw std::invalid_argument("t must be less than or equal to k");
    }

    std::set<int> unique_numbers;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, k - 1);

    // 生成 t 个不同的随机整数
    while (unique_numbers.size() < t) {
        unique_numbers.insert(dis(gen));
    }

    // 将 set 转换为 vector 并返回
    return std::vector<int>(unique_numbers.begin(), unique_numbers.end());
}

