#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double proot(double ct);

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
		                norm_square += x[i]*x[i];
		        }
		        double norm_original = proot(norm_square);

			cout << "norm: " << norm_original << endl;
		}
	}
	fout.close();
}

double proot(double ct){

    double a = ct;
    double b = ct - 1.;
    for(int i=0; i<5; i++){
        a = a * (1. - 0.5 * b);
        b = (b*b) * ((b - 3.)/4.);
    }
    a = a * (1. - 0.5 * b);

    return a;
}
