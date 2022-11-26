#include <iostream>
#include <ostream>
#include <vector>
#include <time.h>
#include <cmath>

double inner_product(std::vector<double> a, std::vector<double> b);

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

	//Calculate the first term
	float time2 = -clock();
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

	time2 += clock();
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "The time of evaluation: " << time2 << " s" << std::endl;

        std::cout << variance << std::endl;

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
