#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
	int n = 1; //The number of data

	ofstream fout;
        fout.open("norm_plaintext.txt"); //Save the results of time
	
        for (int d = 14; d < 30; d++) { //dimension
		std::vector<double> result;
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

	        	//Calculate original norm
        		double norm_square = 0.0;
			std::vector<double> x = dataset_original[0];

		        for (int i = 0; i < d; i++) {
		                norm_square += pow(x[i], 2);
		        }
		        double norm_original = sqrt(norm_square);

			cout << "norm: " << norm_original << endl;
		}
	}
	fout.close();
}


