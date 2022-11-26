#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include <numeric>

using namespace std;

int main() {
	//Setting
	int n = 15; //The number of data

	for (int d = 14; d < 30; d++) {
		std::vector<double> result;
		cout << "d: " << d << " starts" << endl;
		for (int testcase = 0; testcase < 1; testcase++) {
			//Generate dataset
			std::vector<std::vector<double>> dataset(n, std::vector<double> (d, 0.0));

	                for (int z = 0; z < n; z++) {
        	                   std::vector<double> x;
				   double data = (double)z;
                	           int count = 0;
                        	   while (count < d) {
	                                   x.push_back(data);
					   data += 1;
        	                           count += 1;     
				   }     
				   dataset[z] = x;
			}
	
			//Calculate mean vector
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

			cout << "variance: " << variance << endl;
		}
	}
}

