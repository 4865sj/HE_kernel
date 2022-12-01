#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <numeric>

using namespace std;
using namespace chrono;

double inner_product(std::vector<double> a, std::vector<double> b);
double calculator_mean(std::vector<double> input);
double calculator_standard_deviation(std::vector<double> input);
double proot(double ct);

int main() {
	//Setting
	int n = 2; //The number of data
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(10, 20); //Random distribution

	int trials = 20000;

        ofstream fout;
        fout.open("similarity_kernel_plaintext.txt"); //Save the results of time

        for (int d = 4; d < 16; d++) { //dimension
		std::vector<double> result1;
                std::vector<double> result2;
		cout << "d: " << d << " starts with trials " << trials << endl;
                for (int testcase = 0; testcase < trials; testcase++) { //Test trials case for each dimension
                        //Generating dataset
                        std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));

                        for (int z = 0; z < n; z++) {
                                std::vector<double> x;
                                int count = 0;
                                while (count < d) {
                                        double data = dist(gen);
                                        x.push_back(data);
                                        count += 1;
                                }
                                dataset_original[z] = x;
                        }

			//Make kernel
			system_clock::time_point start_time1 = system_clock::now();

			std::vector<std::vector<double>> kernel(n, std::vector<double> (d, 0));
	
			for (int i = 0; i < n; i++) {
				kernel[i][i] = inner_product(dataset_original[i], dataset_original[i]);
			} //Diagonal elements

			for (int i = 0; i < n; i++) {
				for (int j = i+1; j < n; j++) {
					kernel[i][j] = inner_product(dataset_original[i], dataset_original[j]);
					kernel[j][i] = inner_product(dataset_original[i], dataset_original[j]);
				}
			} //Non-diagonal elements

			system_clock::time_point end_time1 = system_clock::now();
			nanoseconds nano1 = end_time1 - start_time1;

		        //Calculate original similarity
		        double denominator = 0.0, numerator = 0.0;

			system_clock::time_point start_time2 = system_clock::now();

			numerator = kernel[0][1];
			denominator = proot(kernel[0][0]*kernel[1][1]);
        
 		        double similarity_original = numerator/denominator;
			
			system_clock::time_point end_time2 = system_clock::now();
			nanoseconds nano2 = end_time2 - start_time2;
			similarity_original = similarity_original;

			result1.push_back(nano1.count());
                        result2.push_back(nano2.count());
		}
		double mean1 = calculator_mean(result1);
                double standard_deviation1 = calculator_standard_deviation(result1);
                cout << "The mean of time of making kernel: " << mean1 << " and the standard deviation of time of making kernel: " << standard_deviation1 << endl;

                double mean2 = calculator_mean(result2);
                double standard_deviation2 = calculator_standard_deviation(result2);
                cout << "The mean of time of evaluation: " << mean2 << " and the standard deviation of time: " << standard_deviation2 << endl;

	}
	fout.close();

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

double proot(double ct){

    double a = ct;
    double b = ct - 1.;
    for(int i=0; i<5; i++){
        a = a * (1. - 0.5 * b);
        b = (b*b) * ((b - 3.)/4.);
    }
    a = a * (1. - 0.5 * b);

    return a;
}
