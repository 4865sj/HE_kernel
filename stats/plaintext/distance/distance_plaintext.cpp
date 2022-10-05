#include <iostream>
#include <ostream>
#include <vector>
#include <time.h>
#include <cmath>

int main() {
	//Setting
	int d = 4; //dimension

        //Inputs
        std::vector<double> x = {0.3, 0.6, 0.5, 0.7};
        std::vector<double> y = {0.7, 0.4, 0.6, 0.3};

        //Calculate original distance
        double distance_square = 0.0;

	float time = -clock(); //Measure time of evaluation

        for (int i = 0; i < d; i++) {
                distance_square += pow(x[i]-y[i], 2);
        }
        double distance_original = sqrt(distance_square);

	time += clock(); //Measure time of evaluation
	time = time/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time << std::endl;

        std::cout << distance_original << std::endl;

}
