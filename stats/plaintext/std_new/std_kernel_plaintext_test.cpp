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

double proot(double ct);

int main() {
	//Setting
	int n = 15;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

	for (int d = 1; d < 16; d++) {//dimension
		cout << "d: " << d << " starts" << endl;

		//Generating dataset
		vector<vector<double>> dataset_original(n, vector<double> (d, 0.));

		for (int z = 0; z < n; z++) {
			vector<double> x;
			int count = 0;
			while (count < d) {
				double data = dist(gen);
				x.push_back(data);
				count += 1;
			}
			dataset_original[z] = x;
		}

		//Making kernel; Scalar multiplication kernel
		vector<vector<double>> kernel_zero(n, vector<double> (n, 0.));
		vector<vector<vector<double>>> kernel_vector(n, kernel_zero); //The vector of kernels; the number of elements is d

		for (int i = 0; i < n; i++) {//Select data_i
			vector<double> data_i = dataset_original[i];
			for (int j = 0; j < n; j++) {//Select data_j
				vector<double> data_j = dataset_original[j];

				for (int k = 0; k < d; k++) {//Select dimension
					kernel_vector[k][i][j] = data_i[k]*data_j[k]; //Update i,j-th element of k-th kernel
				}
			}
		}

		//Calculate std
		system_clock::time_point start_time = system_clock::now();

		//First term
		vector<double> first_term (d, 0.);
		for (int i = 0; i < d; i++) {//Select kernel & dimension
			vector<vector<double>> kernel_i = kernel_vector[i];
			for (int j = 0; j < n; j++) {//Select row & column
				first_term[i] = first_term[i] + kernel_i[j][j];
			}
		}//first_term is the summation of diagonal elements of kernel, now

		for (int i = 0; i < d; i++) {//Select dimension
			first_term[i] = first_term[i]/(double)n; //Divide by n
		}

		//Second term
		vector<double> second_term (d, 0.);
		for (int i = 0; i < d; i ++) {//Select kernel & dimension
			vector<vector<double>> kernel_i = kernel_vector[i];
			for (int j = 0; j < n; j++) {//Select row
				for (int k = 0; k < n; k++) {//Select column
					second_term[i] = second_term[i] + kernel_i[j][k];
				}
			}
		}//second_term is the summation of all elements of kernel, now

		for (int i = 0; i < d; i++) {//Select dimension
			second_term[i] = second_term[i]/(double)(n*n); //Divide by n^2
		}

		//std vector
		vector<double> std (d, 0.); //Initial value
		for (int i = 0; i < d; i++) {//Select dimension
			double variance_i = first_term[i] - second_term[i];
			std[i] = proot(variance_i);
		}

		system_clock::time_point end_time = system_clock::now(); //Measure time of evaluation
		nanoseconds nano = end_time - start_time;

		cout << "The time of evaluation: " << nano.count() << endl;
	}

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






