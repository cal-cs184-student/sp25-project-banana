#include <random>

inline double random_uniform() {
    static std::minstd_rand rng(1234567);
    static std::uniform_real_distribution<double> U(0.0, 1.0);
    return U(rng);
}