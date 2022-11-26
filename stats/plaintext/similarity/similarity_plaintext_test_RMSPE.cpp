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
        fout.open("similarity_plaintext.txt"); //Save the results of time

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

			//Calculate original similarity
		        double dot = 0.0, x_norm_square = 0.0, y_norm_square = 0.0;
			std::vector<double> x = dataset_original[0];
			std::vector<double> y = dataset_original[1];

		        for (int i = 0; i < d; i++) {
                		dot += x[i]*y[i];
		                x_norm_square += x[i]*x[i];
		                y_norm_square += y[i]*y[i];
		        }
		        double similarity_original = dot/(sqrt(x_norm_square)*sqrt(y_norm_square));

			cout << "similarity: " << similarity_original << endl;
		}
	}
	fout.close();

}

