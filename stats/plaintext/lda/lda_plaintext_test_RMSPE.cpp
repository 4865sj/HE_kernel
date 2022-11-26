#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>

using namespace std;
using namespace chrono;

vector<vector<double>> pmatrix_mult(vector<vector<double>> m1,vector<vector<double>> m2, int n, int d, int q);
vector<double> pmatrix_vector_mult(vector<vector<double>> m1,vector<double> m2, int n, int d);
vector<vector<double>> ptranspose(vector<vector<double>> m, int n, int d);
vector<vector<double>> pmatrix_inverse(vector<vector<double>> A);
double power_method(vector<vector<double>> A);

int main() {
	//Setting
        int n = 10; //The number of data
        int c = 3; //The number of classes
        std::vector<int> classes = {3, 3, 4}; //The number of data per class; In this case, The number of data in first class is 3

	for (int d = 4; d < 16; d++) { //Dimension

		std::vector<double> result;
                cout << "d: " << d << " starts" << endl;

		std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0.0));

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


		vector<vector<double>> epsilon(d, vector<double> (d, 0.));
		for (int i = 0; i < d; i++) {
			epsilon[i][i] = 0.001;
		}

		std::vector<std::vector<double>> mean_vector_plaintext(c, std::vector<double> (d, 0.0));

	        int count = 0;

	        for (int i = 0; i < c; i++) { //Select class
              		std::vector<double> mean_plaintext(d, 0.0);
		        for (int j = 0; j < d; j++) {
                	        double tempt = 0.0;
	                        for (int k = 0; k < classes[i]; k++) {
               		                tempt = tempt + dataset_original[count + k][j];
	                        }
               		        tempt = tempt/(double)classes[i];
	                        mean_plaintext[j] = tempt;
               		}
	                mean_vector_plaintext[i] = mean_plaintext;
               		count = count + classes[i];
	        }

	        std::vector<std::vector<double>> A_plaintext(d, std::vector<double> (d, 0.0));

		for (int i = 0; i < c - 1; i++) {
               		for (int j = i + 1; j < c; j++) {
	                        std::vector<double> delta_plaintext(d, 0.0);
               		        for (int k = 0; k < d; k++) {
                               		delta_plaintext[k] = mean_vector_plaintext[i][k] - mean_vector_plaintext[j][k];
	                        } //u_i - u_j

               		        for (int p = 0; p < d; p++) {
                               		A_plaintext[p][p] = A_plaintext[p][p] + delta_plaintext[p]*delta_plaintext[p];
	                        }

               		        for (int p = 0; p < d - 1; p++) {
                               		for (int q = p + 1; q < d; q++) {
	                                        A_plaintext[p][q] = A_plaintext[p][q] + delta_plaintext[p]*delta_plaintext[q];
               		                }
	                        }
               		}
	        }

		for (int p = 0; p < d - 1; p++) {
               		for (int q = p + 1; q < d; q++) {
	                        A_plaintext[q][p] = A_plaintext[p][q];
               		}
	        }

	        std::vector<std::vector<double>> B_plaintext(d, std::vector<double> (d, 0.0));

	        count = 0;

	       	for (int i = 0; i < c; i++) {
			std::vector<std::vector<double>> Z_plaintext(classes[i], std::vector<double> (d, 0.0));
	                for (int p = 0; p < classes[i]; p++) {
               		        for (int q = 0; q < d; q++) {
                               		Z_plaintext[p][q] = dataset_original[count + p][q] - mean_vector_plaintext[i][q];
	                        }
               		}

			std::vector<std::vector<double>> Z_trans_plaintext(d, std::vector<double> (classes[i], 0.0));
			for (int p = 0; p < classes[i]; p++) {
				for (int q = 0; q < d; q++) {
					Z_trans_plaintext[q][p] = Z_plaintext[p][q];
				}
			}

			std::vector<std::vector<double>> B_i_plaintext = pmatrix_mult(Z_trans_plaintext, Z_plaintext, d, classes[i], d);

               		for (int p = 0; p < d; p++) {
				for (int q = 0; q < d; q++) {
					B_plaintext[p][q] = B_plaintext[p][q] + B_i_plaintext[p][q];
				}
			}

			count = count + classes[i];
	        }
			
		for (int i = 0; i < d; i++) {
			B_plaintext[i][i] = B_plaintext[i][i] + epsilon[i][i];
		}

		vector<vector<double>> B_plaintext_inverse = pmatrix_inverse(B_plaintext);
		vector<vector<double>> Final_matrix_plaintext = pmatrix_mult(B_plaintext_inverse, A_plaintext, d, d, d);

		double eigenvalue = power_method(Final_matrix_plaintext);
		cout << "Eigenvalue: " << eigenvalue << endl;
	}
}



vector<vector<double>> pmatrix_mult(vector<vector<double>> m1,vector<vector<double>> m2, int n, int d, int q){

    vector<vector<double>> m;
    vector<double> tmp;
    for(int p=0; p<n; p++){
        for(int i=0; i<q; i++){
            double zero = 0.;
            for(int j=0; j<d; j++){
                zero += m1[p][j] * m2[j][i];
            }
            tmp.push_back(zero);
        }
        m.push_back(tmp);
        tmp.clear();
    }

    return m;
}

vector<double> pmatrix_vector_mult(vector<vector<double>> m1,vector<double> m2, int n, int d){

    vector<double> tmp;
    for(int p=0; p<n; p++){
        double zero = 0.;
        for(int i=0; i<d; i++){
           zero += m1[p][i] * m2[i];
        }
        tmp.push_back(zero);
    }

    return tmp;
}

vector<vector<double>> ptranspose(vector<vector<double>> m, int n, int d){

    vector<vector<double>> transpose(d,vector<double>(n,0.));
    for(int i=0; i<d; i++){
        for(int j=0; j<n; j++){
            transpose[i][j] = m[j][i];
        }
    }

    return transpose;
}

double pinverse(double ct){

    double c = 1. - ct;
    double v = 2. - ct;

    for(int i=0; i<5; i++){
        c = c * c;
        v = v * (1. + c);
    }

    return v;
}

vector<vector<double>> pmatrix_inverse(vector<vector<double>> A) {
	//Input: A is nxn matrix
	
	//Step 1) Setting & calculate C
	int n = A.size();
	vector<vector<double>> A_trans = ptranspose(A, n, n);
	vector<vector<double>> C = pmatrix_mult(A, A_trans, n, n, n);

	//Step 2) Iteration part for getting the dominant eigenvalue
	
	//The first procedure of iteration part; x_new[i] = C[i][0] + C[i][1] + ...
	vector<double> x_new(n, 0.0);
	for (int i = 0; i < n; i++) { //Select row
		double increment = 0.0;
		for (int j = 0; j < n; j++) { //Select column
			increment += C[i][j];
		}
		x_new[i] += increment; //Update x_new[i]
	}

	//Other procedures of interation part
	vector<double> x_old = x_new;
	for (int i = 0; i < 4; i++) { //Iteration with 4 times except first procedure
		x_old = x_new;
		x_new = pmatrix_vector_mult(C, x_old, n, n);
	}

	double alpha_denominator = alpha_denominator; //For eliminating error; this variable is not used
	alpha_denominator = pinverse(x_new[0]); //Because of domain of piverse function, we have some limitation for using this function directly. However, we should evaluate these value for comparing time between for ciphertext and plaintext
	double alpha = x_old[0]/x_new[0]; //alpha = x_old[0]/x_new[0] is the inverse of the dominant eigenvalue

	//Step 3) Iteration part for getting inverse matrix of A
	
	//Generating the necessary matrix for calculating
	vector<vector<double>> alpha_I(n, vector<double> (n, 0.0));
	for (int i = 0; i < n; i++) {
		alpha_I[i][i] = alpha; //All diagonal elements are alpha
	}//alpha_I = (alpha)*I

	vector<vector<double>> two_I(n, vector<double> (n, 0.0));
	for (int i = 0; i < n; i++) {
		two_I[i][i] = 2.0; //All diagonal elements are 2
	}//two_I = 2*I

	//Iteration part
	vector<vector<double>> B = pmatrix_mult(alpha_I, A_trans, n, n, n); //Initial value

	for (int i = 0; i < 15; i++) { //Iteration with 15 times
		//First, calculate latter matrix; 2*I - A*B
		vector<vector<double>> latter_matrix = two_I;
		vector<vector<double>> AB = pmatrix_mult(A, B, n, n, n);
		for (int j = 0; j < n; j++) { //Select row
			for (int k = 0; k < n; k++) { //Select column
				latter_matrix[j][k] = latter_matrix[j][k] - AB[j][k];
			}
		}

		//Update B
		B = pmatrix_mult(B, latter_matrix, n, n, n);
	}

	return B; //The output is the inverse of A
}

double power_method(vector<vector<double>> A) {
        //Input: A is nxn matrix

        //Step 1) Setting
        int n = A.size();

	//Step 2) Iteration part for getting the dominant eigenvalue

        //The first procedure of iteration part; x_new[i] = A[i][0] + A[i][1] + ...
        vector<double> x_new(n, 0.0);
        for (int i = 0; i < n; i++) { //Select row
                double increment = 0.0;
                for (int j = 0; j < n; j++) { //Select column
                        increment += A[i][j];
                }
                x_new[i] += increment; //Update x_new[i]
        }

        //Other procedures of interation part
        vector<double> x_old = x_new;
        for (int i = 0; i < 4; i++) { //Iteration with 4 times except first procedure
                x_old = x_new;
                x_new = pmatrix_vector_mult(A, x_old, n, n);
        }

	double lambda_denominator = lambda_denominator; //For eliminate error; this variable is not used
	lambda_denominator = pinverse(x_old[0]); //Because of domain of piverse function, we have some limitation for using this function directly. However, we should evaluate these value for comparing time between for ciphertext and plaintext
	double lambda = x_new[0]/x_old[0]; //lambda = x_new[0]/x_old[0] is the dominant eigenvalue

	return lambda; //The output is the dominant eigenvalue of A
}


