#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
	//Setting
	int n = 2; //The number of data

	ofstream fout;
        fout.open("distance_plaintext.txt"); //Save the results of time

        for (int d = 14; d < 30; d++) { //dimension
                cout << "d: " << d << " starts" << endl;
                for (int testcase = 0; testcase < 1; testcase++) { //Test 1 case for each dimension
                        //Generating dataset
                        std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));

                        for (int z = 0; z < n; z++) {
                                   std::vector<double> x;
                                   double data = (double)z;
                                   int count = 0;
                                   while (count < d) {
                                           x.push_back(data);
                                           data += 1;
                                           count += 1;
                                   }
                                   dataset_original[z] = x;
                        }

			//Calculate original distance
		        double distance_square = 0.0;
			std::vector<double> x = dataset_original[0];
			std::vector<double> y = dataset_original[1];

		        for (int i = 0; i < d; i++) {
		                distance_square += pow(x[i]-y[i], 2);
		        }
		        double distance_original = sqrt(distance_square);

			cout << "distance: " << distance_original << endl;
		}
	}
	fout.close();
}

