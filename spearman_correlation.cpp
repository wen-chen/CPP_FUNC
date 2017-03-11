#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

std::vector <double> GetRanks(std::vector <double> & data) {
	int len = data.size();
	std::vector <int> index(len, 0);
	std::vector <double> rank(len, 0);
	for (int i = 0; i < len; ++i) {
		index[i] = i;
	}
	std::sort(index.begin(), index.end(),
	    [&](const int & a, const int & b) {
	    	return (data[a] > data[b]);
		}
	);

    int sumranks = 0;
    int dupcount = 0;
    for (int i = 0; i < len; ++i) {
    	sumranks = sumranks + i;
    	dupcount = dupcount + 1;
    	if ((i == len -1) || (data[index[i]] != data[index[i + 1]]) ) {
    		double averank = double(sumranks) / double(dupcount) + 1;
    		for (int j = i - dupcount + 1; j < i + 1; ++j) {
    			rank[index[j]] = averank;
			}
			sumranks = 0;
			dupcount = 0;
		}
	}
	return rank;
}

double pearson (std::vector <double> & v1, std::vector <double> & v2) {
	int n = v1.size();
	
	double mean_1 = 0;
	double mean_2 = 0;
	for (int i=0; i < n; i++) {
		mean_1 += v1[i];
		mean_2 += v2[i];
	}		
	mean_1 = mean_1 / n;
	mean_2 = mean_2 / n;

	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	for (int i=0; i < n; i++) {
		v1[i] = v1[i] - mean_1;
		v2[i] = v2[i] - mean_2;
		sum1 += v1[i] * v2[i];
		sum2 += v1[i] * v1[i];
		sum3 += v2[i] * v2[i];
	}
	
	double cor = sum1 / std::sqrt(sum2 * sum3);
	return cor;	
}

double spearman(std::vector <double> & v1, std::vector <double> & v2) {
    std::vector <double> v1_rank = GetRanks(v1);
    std::vector <double> v2_rank = GetRanks(v2);
    return pearson(v1, v2);
}

const double STOP = 1.0e-8;
const double TINY = 1.0e-30;
double incbeta(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 1.0/0.0;

    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a + 1.0) / (a + b + 2.0)) {
        return (1.0 - incbeta(b, a, 1.0 -x )); /*Use the fact that beta is symmetrical.*/
    }

    /*Find the first part before the continued fraction.*/
    const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
    const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    double f = 1.0, c = 1.0, d = 0.0;

    int i, m;
    for (i = 0; i <= 200; ++i) {
        m = i/2;

        double numerator;
        if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        } else if (i % 2 == 0) {
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
        } else {
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
        }

        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (fabs(d) < TINY) d = TINY;
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (fabs(c) < TINY) c = TINY;

        const double cd = c*d;
        f *= cd;

        /*Check for stop.*/
        if (fabs(1.0-cd) < STOP) {
            return front * (f-1.0);
        }
    }

    return 1.0/0.0; /*Needed more loops, did not converge.*/
}

double student_t_cdf(const double t, const double v) {
    double x = (t + std::sqrt(t * t + v)) / (2.0 * std::sqrt(t * t + v));
    double prob = incbeta(v * 0.5, v * 0.5, x);
    return prob;
}

double get_p_value(double r, double n) {
	double t = std::abs(r) * std::sqrt((n -2.0) / (1 - r * r));
	double p = (1.0 - student_t_cdf(t, n - 2.0)) * 2.0;
	return p;	
}

int main() {
    std::vector <double> vec1 = {56.0,75.0,45.0,71.0,61.0,64.0,58.0,80.0,76.0,61.0};
	std::vector <double> vec2 = {66.0,70.0,40.0,60.0,65.0,56.0,59.0,77.0,67.0,63.0};
	double corr = spearman(vec1, vec2);
	double p = get_p_value(corr, 5);
    std::cout << corr << std::endl;
	std::cout << p << std::endl;
    return 0;
}