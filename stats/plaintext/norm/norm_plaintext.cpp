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

        //Calculate original norm
        double norm_square = 0.0;

	float time = -clock(); //Measure time of evaluation

        for (int i = 0; i < d; i++) {
                norm_square += pow(x[i], 2);
        }
        double norm_original = sqrt(norm_square);

	time += clock(); //Measure time of evaluation
	time = time/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time << std::endl;

        std::cout << norm_original << std::endl;

}
