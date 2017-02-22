#include <cmath>
#include <iostream>

const double PI = 3.141592653589793238463;

double t_pdf(double n, double t) {
	return std::tgamma((n + 1) * 0.5) * std::pow(1 + t * t / n, (n + 1) * -0.5) / 
    (std::sqrt(n * PI) * std::tgamma(n * 0.5));
}

int main() {
	double t = t_pdf(8.0, 0);
	std::cout << t << std::endl;
	return 0;
}


