#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

using namespace std;
using namespace chrono;

int main() {
	int n = 1; //The number of data
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(10, 20); //Random distribution

        ofstream fout;
        fout.open("norm_plaintext.txt"); //Save the results of time
	
        for (int d = 4; d < 401; d = d + 4) { //dimension
        	for (int testcase = 0; testcase < 10; testcase++) { //Test 10 case for each dimension
                        cout << "Test dimension " << d << " and case " << testcase + 1 << endl;
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


	        	//Calculate original norm
        		double norm_square = 0.0;
			std::vector<double> x = dataset_original[0];

			system_clock::time_point start_time = system_clock::now();

		        for (int i = 0; i < d; i++) {
		                norm_square += pow(x[i], 2);
		        }
		        double norm_original = sqrt(norm_square);

			system_clock::time_point end_time = system_clock::now(); //Measure time of evaluation
			nanoseconds nano = end_time - start_time;

			norm_original = norm_original;

			fout << nano.count() << endl;
		}
	}
	fout.close();
}
