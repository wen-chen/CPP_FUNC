#include <iostream>

unsigned long long comb(unsigned long long n, unsigned long long k) {
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {        
        r = r * n;
        r = r / d;
        n = n - 1;
    }
    return r;
}

int main() {
    std::cout << comb(20, 5) << std::endl;
    return 0;
}

