#include "openfhe.h"
#include <time.h>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace lbcrypto;
using namespace std;
using namespace chrono;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
add_many(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, CryptoContext<DCRTPoly> cc);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
transpose(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
diagonalize(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
matrix_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, int n, CryptoContext<DCRTPoly> cc);

vector<vector<double>> pmatrix_mult(vector<vector<double>> m1,vector<vector<double>> m2, int n, int d, int q);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
power_method(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> c, int n, CryptoContext<DCRTPoly> cc);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >>
matrix_inverse(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, double dv, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
matrix_vector_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxt, int n, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc);

int main() {
	ofstream fout;
	fout.open("lda.txt");

	for (int d = 5; d < 16; d++) {
		cout << "d: " << d << " starts" << endl;
	        // Setting
        	uint32_t multDepth = 65;
	        uint32_t scaleModSize = 50;
        	uint32_t batchSize = 64;
		int n = 10; //The number of data
		int c = 3; //The number of classes
		int nd = max(n, d);
		std::vector<int> classes = {3, 3, 4}; //The number of data per class; In this case, The number of data in first class is 2
	        
		CCParams<CryptoContextCKKSRNS> parameters;
        	parameters.SetMultiplicativeDepth(multDepth);
	        parameters.SetScalingModSize(scaleModSize);
        	parameters.SetBatchSize(batchSize);

	        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        	cc -> Enable(PKE);
	        cc -> Enable(KEYSWITCH);
        	cc -> Enable(LEVELEDSHE);

		auto keys = cc -> KeyGen();

        	vector<int> ro;
	        for (int i = 1; i < nd + 1; i++) {
        	        ro.push_back(i);
                	ro.push_back(-i);
	        }

        	cc -> EvalMultKeyGen(keys.secretKey);
	        cc -> EvalRotateKeyGen(keys.secretKey, ro);
	
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


		std::vector<Plaintext> ptxt_one;

        	for (int i = 0; i < d; i++) {
                	std::vector<double> one(d, 0.0);
	                one[i] = 1.0;
        	        Plaintext one_tmp = cc -> MakeCKKSPackedPlaintext(one);
                	ptxt_one.push_back(one_tmp);
	        }

		vector<double> zero(d, 0.0);
		Plaintext ptxt_zero = cc -> MakeCKKSPackedPlaintext(zero); //Encoding zero vector
		auto ct_zero = cc -> Encrypt(keys.publicKey, ptxt_zero); //Encrypting zero vector

		std::vector<Plaintext> plaintext;

		for (int i = 0; i < n; i++) {
			Plaintext ptxt = cc -> MakeCKKSPackedPlaintext(dataset_original[i]);
			plaintext.push_back(ptxt);
			std::cout << "Input x" << i << ": " << ptxt << std::endl;
		}

	        std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dataset;
	
		for (int i = 0; i < n; i++) {
			auto c = cc -> Encrypt(keys.publicKey, plaintext[i]);
			dataset.push_back(c);
		}

                //Generate epsilon
                std::vector<Plaintext> ptxt_epsilon;
                std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > epsilon;
                for (int i = 0; i < n; i++) {
                        std::vector<double> epsilon_row(n, 0.0);
                        epsilon_row[i] = 0.001;
                        Plaintext epsilon_tmp = cc -> MakeCKKSPackedPlaintext(epsilon_row);
                        ptxt_epsilon.push_back(epsilon_tmp);
                        auto ct_epsilon = cc -> Encrypt(keys.publicKey, epsilon_tmp);
                        epsilon.push_back(ct_epsilon);
                }


		//Calculate mean vectors per classes
	
		std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> mean_vector;

	        int count = 0; //Save the number of data in previous classes
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > mean;

		for (int i = 0; i < c; i++) {
			mean = dataset[count];
			int number = classes[i]; //The number of data in class[i]
			for (int j = 1; j < number; j++) {
				mean = cc -> EvalAdd(mean, dataset[count + j]);
			}
			double a = 1.0/(double)number;
        	        mean = cc -> EvalMult(mean, a); //mean is class[i]'s mean vector
                	mean_vector.push_back(mean);//mean_vector[i] is class[i]'s mean vector
			count = count + classes[i];
		}
		cout << "mean vector complete" << endl;
	
		//Calculate A

		std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > A(d, ct_zero); //ct_zero is just temparary data
	
		for (int i = 0; i < c - 1; i++) { //Select u_i
                	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > u_i = mean_vector[i];
	                for (int j = i + 1; j < c; j++) { //Select u_j
        	                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > u_j = mean_vector[j];
                	        //First, calculate diagonal elements
                        	for (int k = 0; k < d; k++) { //Select row
                                	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta = cc -> EvalSub(u_i, u_j);
	                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment = cc -> EvalMult(delta, delta); //kth element is increment of M[k][k]
        	                        increment = cc -> EvalMult(increment, ptxt_one[k]); //Reset data except kth element
                	                A[k] = cc -> EvalAdd(A[k], increment); //Update M[k][k]
                        	} //Diagonal elements
	                        //Next, calculate non-diagonal elements
        	                for (int k = 0; k < d - 1; k++) { //Select row
                	                for (int l = k + 1; l < d; l++) { //Select column
                        	                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta1 = cc -> EvalSub(u_i, u_j); //We need kth element
                                	        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta1_rotate = cc -> EvalRotate(delta1, -(l - k)); //delta1_rotate[l] = delta1[k]

                                        	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta2 = cc -> EvalSub(u_i, u_j); //We need lth element

	                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment = cc -> EvalMult(delta1_rotate, delta2); //lth element is increment of A[k][l]; delta[k]*delta[l]
        	                                increment = cc -> EvalMult(increment, ptxt_one[l]); //Reset data except lth element
                	                        A[k] = cc -> EvalAdd(A[k], increment); //Update A[k][l]

                        	                increment = cc -> EvalRotate(increment, l - k); //kth element is increment of A[l][k]
                                	        A[l] = cc -> EvalAdd(A[l], increment); // By symmetry, update A[l][k]
	                                }
        	                } //Non-diagonal elements
                	}
	        } //The output is A
		cout << "A complete" << endl;
	
		//Calculate B
	
		std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > B(d, ct_zero); //ct_zero is just temparary data

		count = 0; //Save the number of data in previous classes

		for (int i = 0; i < c; i++) { //Select class
			std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > Z(classes[i], ct_zero); //ct_zero is just temparary data
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment;

			for (int p = 0; p < classes[i]; p++) { //Select data in class
				Z[p] = cc -> EvalSub(dataset[count + p], mean_vector[i]);
			} //Z[p] = d_p - u_i


			int max_d = max(classes[i], d);
			if (classes[i] < max_d) {
				for (int j = 0; j < max_d - classes[i]; j++) {
					Z.push_back(ct_zero);
				}
			} //For using Ahn's functions, we should set the size of Z

			std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > Z_diag = diagonalize(Z, ptxt_one, max_d, cc);
			std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > Z_trans = transpose(Z, ptxt_one, max_d, cc);
			std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > B_i = matrix_mult(Z_diag, Z_trans, max_d, cc);
			for (int p = 0; p < d; p++) {
				B[p] = cc -> EvalAdd(B[p], B_i[p]);
			}		
			count = count + classes[i];
		} //The output is B
		cout << "B complete" << endl;

		//Calculate dominant eignevector of (B_inverse)*(A)
                for (int i = 0; i < d; i++) {
                        B[i] = cc -> EvalAdd(B[i], epsilon[i]);
                } //Since B is singular, we need to add epsilon

		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> B_inverse = matrix_inverse(B, d, 30.0, cc); //The transpose of B_inverse
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> B_inverse_diag = diagonalize(B_inverse, ptxt_one, d, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> A_trans = transpose(A, ptxt_one, d, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> Final_matrix_trans = matrix_mult(B_inverse_diag, A_trans, d, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> Final_matrix = transpose(Final_matrix_trans, ptxt_one, d, cc);
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > eigenvector = power_method(Final_matrix, d, cc);

		//Decrypting and Decoding
		Plaintext result1;

		cc -> Decrypt(keys.secretKey, eigenvector, &result1);
		result1 -> SetLength(d);

		std::cout << result1 << std::endl; //Print eigenvector over encrypted data
		
	}
	fout.close();
			
}

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
transpose(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc){

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > ctxst;
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            auto tmp = cc->EvalMult(ctxs[j],diag[i]);
            if(j < i){
                tmp = cc->EvalRotate(tmp,i-j);
            }
            else if(j > i){
                tmp = cc->EvalRotate(tmp,-j+i);
            }
            tmpv.push_back(tmp);
        }
        auto tmp = add_many(tmpv,n,cc);
        tmpv.clear();
        ctxst.push_back(tmp);
    }

    return ctxst;
}

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
diagonalize(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc){
//Input data are ctxst: n x m matrix s.t n >= m, diag: {{1, 0, 0, ...}, {0, 1, 0, ...}, ... , {0, 0, ... , 1}}, n: the number of rows of ctxst, cc: key
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > dev; //Result matrix
    //First, extract the diagonal elements
    for(int i=0; i<n; i++){
        auto tmp = cc->EvalMult(ctxst[i],diag[i]); //Extract nth diagonal element
        tmpv.push_back(tmp);
    }
    auto tmp = add_many(tmpv,n,cc); //Combine extractions
    dev.push_back(tmp);
    tmpv.clear();

    //Next, extract the non-diagonal elements
    for(int i=1; i<n; i++){ //Select the row of result matrix
        for(int j=0; j<n; j++){ //Select the row of input matrix
            auto tmp = cc->EvalMult(ctxst[j],diag[(n+j-i)%n]); //dev[i][(n+j-i)%n] = ctxst[j][(n+j-i)%n]
            tmpv.push_back(tmp);
        }
        auto tmp = add_many(tmpv,n,cc); //Combine extractions, for each i
        dev.push_back(tmp);
        tmpv.clear();
    }

    return dev;
} //The output is the transpose of ctxst's diag matrix; 1st row is diagonal elements of input matrix, 2nd row is -1 shift diagonal elements of input matrix, ...

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
matrix_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, int n, CryptoContext<DCRTPoly> cc){
//Suppose that we want to calculate a multiplication between p x q matrix and q x r matrix; Call A and B, repectively
//For using this function, we should set the number of row of A and B to max(p, q, r), respectively
//For example, if p is the biggest number, we should transform the second matrix to p x r matrix, by adding zero rows
//Finally, input data are dev: A's diag matrix, ctxst: the transpose of B, n: max(p, q, r), cc: key
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > matrix_mult; //Result matrix
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(j==0){
                auto tmp = cc->EvalMult(dev[j],ctxst[i]);
                tmpv.push_back(tmp);
            }
            else{
                auto tmp = cc->EvalRotate(ctxst[i],j);
                auto tmp1 = cc->EvalRotate(ctxst[i],-n+j);
                tmp1 = cc->EvalAdd(tmp,tmp1);
                tmp = cc->EvalMult(dev[j],tmp1);
                tmpv.push_back(tmp);
            }
        }
        auto tmp = add_many(tmpv,n,cc);
        tmpv.clear();
        matrix_mult.push_back(tmp);
    }

    return matrix_mult; //The output is the transpose of A*B
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
add_many(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, CryptoContext<DCRTPoly> cc) {

    auto res = cc->EvalAdd(ctxs[0],ctxs[1]);
    if (n > 2){
        for(int i=2; i<n; i++){
            res = cc->EvalAdd(res,ctxs[i]);
        }
    }
    return res;
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

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >>
matrix_inverse(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, double dv, CryptoContext<DCRTPoly> cc){
//Suppose that we input nxn matrix A. Then, the output of this function is the transpose of the inverse of A
//We should input some proper divisor, because this function use inverse function, and the domain of inverse function is 0 < x < 1
//Finally, input data are ctxs: nxn matrix, n: row & col, dv: proper divisor, cc: key
    vector<Plaintext> p;
    for(int i=0; i<n; i++){
        Plaintext ptxt;
        p.push_back(ptxt);
    }
    vector<Plaintext> diag;
    for(int i=0; i<n; i++){
        vector<double> d;
        for(int j=0; j<n; j++){
            if(i == j){
                d.push_back(1.);
            }
            else{
                d.push_back(0.);
            }
        }
        p[i] = cc->MakeCKKSPackedPlaintext(d);
        diag.push_back(p[i]);
        d.clear();
    }//diag is the identity matrix

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;

    auto ctxst = transpose(ctxs,diag,n,cc);//The transpose of A

    auto c = matrix_mult(diagonalize(ctxst,diag,n,cc),ctxs,n,cc); //Multiplication between A and the transpose of A
    auto c_diag = diagonalize(c,diag,n,cc); 

    auto tmp = add_many(c,n,cc);
    tmpv.push_back(tmp); //The first procedure of iteration
    for(int i=0; i<4; i++){
        auto tmp = matrix_vector_mult(c_diag,tmpv[i],n,cc); //x_new = C*x_old
        tmpv.push_back(tmp);
    }

    auto binverse = cc->EvalMult(tmpv[4],1./dv); //Rescaling for using inverse function
    binverse = cc->EvalMult(binverse,diag[0]); //Eliminate except the first element
    auto alpha = cc->EvalMult(inverse(binverse,cc),tmpv[3]); //x_new(1)/x_old(1)
    alpha = cc->EvalMult(alpha,1./dv); //Return to original value from rescaling value
    alpha = cc->EvalMult(alpha,diag[0]); //Eliminate except the first element
    tmpv.clear();
    tmpv.push_back(alpha); //The inverse of the dominant eigenvalue
    for(int i=1; i<n; i++){
        auto tmp = cc->EvalRotate(alpha,-i);
        tmpv.push_back(tmp);
    }
    alpha = add_many(tmpv,n,cc);

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > b;
    for(int i=0; i<n; i++){
        auto tmp = cc->EvalMult(alpha,ctxs[i]);
        b.push_back(tmp);
    }//b = (alpha)*(the transpose of A)

    vector<Plaintext> diag2;
    for(int i=0; i<n; i++){
        vector<double> d;
        for(int j=0; j<n; j++){
            if(i == j){
                d.push_back(2.);
            }
            else{
                d.push_back(0.);
            }
        }
        p[i] = cc->MakeCKKSPackedPlaintext(d);
        diag2.push_back(p[i]);
        d.clear();
    }//diag2 is 2*I

    for(int i=0; i<15; i++){
        auto ab = matrix_mult(diagonalize(ctxst,diag,n,cc),b,n,cc); //A*B
        for(int j=0; j<n; j++){
            ab[j] = cc->EvalSub(diag2[j],ab[j]); //2*I - A*B
        }
        b = matrix_mult(diagonalize(b,diag,n,cc),ab,n,cc); //The transpose of B*(2*I - A*B)
    }

    return b; //The output is the transpose of the inverse of A
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
power_method(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> c, int n, CryptoContext<DCRTPoly> cc){

    vector<Plaintext> p;
    for(int i=0; i<n; i++){
        Plaintext ptxt;
        p.push_back(ptxt);
    }
    vector<Plaintext> diag;
    for(int i=0; i<n; i++){
        vector<double> d;
        for(int j=0; j<n; j++){
            if(i == j){
                d.push_back(1.);
            }
            else{
                d.push_back(0.);
            }
        }
        p[i] = cc->MakeCKKSPackedPlaintext(d);
        diag.push_back(p[i]);
        d.clear();
    }

    auto c_diag = diagonalize(c,diag,n,cc);

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;

    auto tmp = add_many(c,n,cc);
    tmpv.push_back(tmp);
    for(int i=0; i<4; i++){
        auto tmp = matrix_vector_mult(c_diag,tmpv[i],n,cc);
        tmpv.push_back(tmp);
    }

    return tmpv[4];
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
matrix_vector_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxt, int n, CryptoContext<DCRTPoly> cc){

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;
    for(int j=0; j<n; j++){
        if(j==0){
            auto tmp = cc->EvalMult(dev[j],ctxt);
            tmpv.push_back(tmp);
        }
        else{
            auto tmp = cc->EvalRotate(ctxt,j);
            auto tmp1 = cc->EvalRotate(ctxt,-n+j);
            tmp1 = cc->EvalAdd(tmp,tmp1);
            tmp = cc->EvalMult(dev[j],tmp1);
            tmpv.push_back(tmp);
        }
    }
    auto matrix_vector_mult = add_many(tmpv,n,cc);
    tmpv.clear();

    return matrix_vector_mult;

}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc) {

    auto ct1 = cc->EvalSub(1.0,ct);
    auto ct2 = cc->EvalAdd(ct1,1.0);

    for(int i=0; i<5; i++){
        ct1 = cc->EvalMult(ct1,ct1);
        ct2 = cc->EvalMult(ct2,cc->EvalAdd(ct1,1.0));
    }

    return ct2;
}
