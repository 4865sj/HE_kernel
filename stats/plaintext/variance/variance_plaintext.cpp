#include <iostream>
#include <ostream>
#include <vector>
#include <time.h>
#include <cmath>

int main() {
	//Setting
	int d = 4; //dimension
	int n = 4; //The number of data

        //Inputs
        std::vector<double> x1 = {2.0, 2.1, 2.2, 2.3};
        std::vector<double> x2 = {3.1, 3.2, 3.3, 3.4};
	std::vector<double> x3 = {4.1, 4.2, 4.3, 4.4};
	std::vector<double> x4 = {5.1, 5.2, 5.3, 5.4};
	std::vector<std::vector<double>> dataset(n, std::vector<double> (d, 0));
	dataset = {x1, x2, x3, x4};

	//Calculate mean vector
	float time = -clock(); //Measure time of evaluation

	std::vector<double> mean;

	for (int i = 0; i < d; i++) {
		mean.push_back(0.0);
	}

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
			variance_square += pow(mean[j] - dataset[i][j], 2);
		}
	}

	double variance = variance_square/(double)n;

	time += clock(); //Measure time of evaluation
	time = time/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time << std::endl;

        std::cout << variance << std::endl;

}

