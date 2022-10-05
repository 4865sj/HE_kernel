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

        //Calculate original similarity
        double dot = 0.0, x_norm_square = 0.0, y_norm_square = 0.0;

	float time = -clock(); //Measure time of evaluation

        for (int i = 0; i < d; i++) {
                dot += x[i]*y[i];
                x_norm_square += x[i]*x[i];
                y_norm_square += y[i]*y[i];
        }
        double similarity_original = dot/(sqrt(x_norm_square)*sqrt(y_norm_square));

	time += clock(); //Measure time of evaluation
	time = time/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time << std::endl;

        std::cout << similarity_original << std::endl;

}
