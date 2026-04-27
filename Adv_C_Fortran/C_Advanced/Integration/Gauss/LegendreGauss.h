// C++ class generated from PHP code
#ifndef LEGENDRE_GAUSS_H
#define LEGENDRE_GAUSS_H

#include <vector>
#include <map>

class LegendreGauss {
public:
    static const std::map<int, std::vector<double>>& get_roots();
    static const std::vector<double>& get_roots(int n);
    static void initialize();

private:
    static std::map<int, std::vector<double>> legendre_roots;
    static bool initialized;
};
#endif // LEGENDRE_GAUSS_H