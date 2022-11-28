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

double inner_product(std::vector<double> a, std::vector<double> b);
double calculator_mean(std::vector<double> input);
double calculator_standard_deviation(std::vector<double> input);

int main() {
        //Setting
        int n = 15; //The number of data

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(10, 20); //Random distribution

	int trials = 20000;

	for (int d = 14; d < 29; d++) {
		std::vector<double> result1;
		std::vector<double> result2;
                cout << "d: " << d << " starts with trials " << trials << endl;
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

			//Make kernel
			system_clock::time_point start_time_kernel = system_clock::now();
			std::vector<std::vector<double>> kernel(n, std::vector<double> (n, 0));
	
			for (int i = 0; i < n; i++) {
				kernel[i][i] = inner_product(dataset[i], dataset[i]);
			} //Diagonal elements

			for (int i = 0; i < n; i++) {
				for (int j = i+1; j < n; j++) {
					kernel[i][j] = inner_product(dataset[i], dataset[j]);
					kernel[j][i] = inner_product(dataset[i], dataset[j]);
				}
			} //Non-diagonal elements
		
			system_clock::time_point end_time_kernel = system_clock::now();
        	        nanoseconds nano1 = end_time_kernel - start_time_kernel;

			//Calculate the first term
			system_clock::time_point start_time_eval = system_clock::now();

			double first_term = 0.0;
			for (int i = 0; i < n; i++) {
				first_term += kernel[i][i];
			}
			first_term = first_term/(double)n;

			//Calculate the second term
			double second_term = 0.0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					second_term += kernel[i][j];
				}
			}
			second_term = second_term/(double)pow(n, 2);

			//Calculate variance
			double variance = first_term - second_term;

			system_clock::time_point end_time_eval = system_clock::now();
	                nanoseconds nano2 = end_time_eval - start_time_eval;
	        	
			result1.push_back(nano1.count());
			result2.push_back(nano2.count());

			variance = variance;
		}
		double mean1 = calculator_mean(result1);
                double standard_deviation1 = calculator_standard_deviation(result1);
                cout << "The mean of time of making kernel: " << mean1 << " and the standard deviation of time of making kernel: " << standard_deviation1 << endl;

		double mean2 = calculator_mean(result2);
                double standard_deviation2 = calculator_standard_deviation(result2);
                cout << "The mean of time of evaluation: " << mean2 << " and the standard deviation of time: " << standard_deviation2 << endl;

	}
}

double inner_product(std::vector<double> a, std::vector<double> b) {
	int n = a.size();
	int m = b.size();
	double result = 0.0;

	if (n != m) {
		std::cout << "Two input vectors have differenct size!" << std::endl;

		return 0.0;
	}
	else {
		for (int i = 0; i < n; i++) {
			result += a[i]*b[i];
		}

		return result;
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

