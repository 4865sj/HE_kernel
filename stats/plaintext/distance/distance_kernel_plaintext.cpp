#include <iostream>
#include <ostream>
#include <vector>
#include <time.h>
#include <cmath>

double inner_product(std::vector<double> a, std::vector<double> b);

int main() {
	//Setting
	int d = 4; //dimension
	int n = 2; //The number of data

        //Inputs
        std::vector<double> x = {0.3, 0.6, 0.5, 0.7};
        std::vector<double> y = {0.7, 0.4, 0.6, 0.3};
	std::vector<std::vector<double>> dataset(n, std::vector<double> (d, 0));
	dataset = {x, y};

	//Make kernel
	float time1 = -clock();
	std::vector<std::vector<double>> kernel(n, std::vector<double> (d, 0));
	
	for (int i = 0; i < n; i++) {
		kernel[i][i] = inner_product(dataset[i], dataset[i]);
	} //Diagonal elements

	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			kernel[i][j] = inner_product(dataset[i], dataset[j]);
			kernel[j][i] = inner_product(dataset[i], dataset[j]);
		}
	} //Non-diagonal elements

	time1 += clock();
	time1 = time1/CLOCKS_PER_SEC;
	std::cout << "The time of making kernel: " << time1 << " s" << std::endl;

        //Calculate original distance
	float time2 = -clock(); //Measure time of evaluation

	double distance_square = kernel[0][0] + kernel[1][1] - 2*kernel[0][1];

        double distance_original = sqrt(distance_square);

	time2 += clock(); //Measure time of evaluation
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time2 << " s" << std::endl;

        std::cout << distance_original << std::endl;

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
