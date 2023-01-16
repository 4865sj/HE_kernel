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
	std::uniform_real_distribution<> dist(10, 20); //Random distribution

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

		//Calculate mean
		system_clock::time_point start_time = system_clock::now();

		vector<double> mean_vector (d, 0.); //Initial value
		for (int i = 0; i < n; i++) {//Select data_i
			vector<double> data_i = dataset_original[i];
			for (int j = 0; j < d; j++) {//Select dimension
				mean_vector[j] = mean_vector[j] + dataset_original[i][j];
			}
		} //mean_vector is the summation of all data, now

		for (int i = 0; i < d; i ++) {//Select dimension
			mean_vector[i] = mean_vector[i]/(double)n; //Divide by n
		} //mean_vector is the mean vector, now

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






