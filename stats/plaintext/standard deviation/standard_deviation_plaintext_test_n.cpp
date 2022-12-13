#include <iostream>
#include <ostream>
#include <vector>
#include <time.h>
#include <cmath>
#include <chrono>
#include <ctime>
#include <random>
#include <numeric>

using namespace std;
using namespace chrono;

double calculator_mean(std::vector<double> input);
double calculator_standard_deviation(std::vector<double> input);
double proot(double ct);

int main() {
	//Setting
	int d = 1; //Dimension

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(10, 20); //Random distribution

	int trials = 20000;

	for (int n = 1; n < 16; n++) {
		std::vector<double> result;
		cout << "n: " << n << " starts with trials " << trials << endl;
		for (int testcase = 0; testcase < trials; testcase++) {
			std::vector<std::vector<double>> dataset(n, std::vector<double> (d, 0.0));

	                for (int z = 0; z < n; z++) {
        	                   std::vector<double> x;
                	           int count = 0;
                        	   while (count < d) {
		    			   double data = dist(gen);
	                                   x.push_back(data);
        	                           count += 1;     
				   }     
				   dataset[z] = x;
			}
	
			//Calculate mean vector
			system_clock::time_point start_time_eval = system_clock::now();

			std::vector<double> mean(d, 0.);

			for (int i = 0; i < d; i++) {
				double tmp = 0.0;
				for (int j = 0; j < n; j++) {
					tmp += dataset[j][i];
				}
				mean[i] = tmp/(double)n;
			}

			//Calculate variance
	
			double variance_square = 0.0;
	
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < d; j++) {
					variance_square += (mean[j] - dataset[i][j])*(mean[j] - dataset[i][j]);
				}
			}

			double variance = variance_square/(double)n;

			//Calculate standard deviation	
			double standard_deviation = proot(variance);

			system_clock::time_point end_time_eval = system_clock::now();
			nanoseconds nano = end_time_eval - start_time_eval;

			result.push_back(nano.count());

			standard_deviation = standard_deviation;
		}
		double mean = calculator_mean(result);
		double standard_deviation = calculator_standard_deviation(result);
		cout << "The mean of time of evaluation: " << mean << " and the standard deviation of time: " << standard_deviation << endl;
	}
}

double calculator_mean(std::vector<double> input) {
	int n = input.size();
	double sum = accumulate(input.begin(), input.end(), 0.);
	double mean = sum/double(n);
	return mean;
}

double calculator_standard_deviation(std::vector<double> input) {
	int n = input.size();
	double mean = calculator_mean(input);
	double sum = 0.;

	for (int i = 0; i < n; i++) {
		double increment = pow(input[i] - mean, 2);
		sum = sum + increment;
	}

	double standard_deviation = sqrt(sum/double(n));

	return standard_deviation;
}

double proot(double ct){

    double a = ct;
    double b = ct - 1.;
    for(int i=0; i<6; i++){
        a = a * (1. - 0.5 * b);
        b = (b*b) * ((b - 3.)/4.);
    }
    a = a * (1. - 0.5 * b);

    return a;
}

