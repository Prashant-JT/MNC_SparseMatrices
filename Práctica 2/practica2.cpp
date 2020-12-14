#include <cstdio>
#include <cstdlib>
#include <string>

#include <mkl.h>

using namespace std;

class Utils {
public:
	void printMatrix(double* mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (mat[i] < 0) {
				printf("%4.1f\t", mat[i]);
			}
			else {
				printf("%4.2f\t", mat[i]);
			}
			if ((i + 1) % n == 0) {
				printf("\n\t\t");
			}
		}
		printf("\n");
	}

	void printVectorInt(int * mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (i == (m * n) - 1) {
				printf("%d\n", mat[i]);
			} else {
				printf("%d, ", mat[i]);
			}			
		}
	}

	void printVectorDouble(double* mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (i == (m * n) - 1) {
				printf("%4.2f\n", mat[i]);
			}
			else {
				if (mat[i] < 0) {
					printf("%4.1f, ", mat[i]);
				}
				else {
					printf("%4.2f, ", mat[i]);
				}
			}
		}
	}

	int countNNZ(double * mat, int size) {
		int count = 0;
		for (int i = 0; i < size; i++) {
			if (mat[i] != 0) {
				count++;
			}
		}
		return count;
	}
	
	double* clone(double* matrix, int size) {
		double* res = new double[size];
		memcpy(res, matrix, sizeof(double) * size);
		return res;
	}

	int* cloneInt(int* matrix, int size) {
		int* res = new int[size];
		memcpy(res, matrix, sizeof(int) * size);
		return res;
	}

	double* processMatrix(string mat, int size) {
		double* vect = new double[size];
		int p = 0;
		string temp = "";

		for (int i = 0; i < mat.size(); i++) {
			if (mat[i] != ' ') {
				temp += mat[i];
			}
			else {
				vect[p] = stof(temp);
				p++;
				temp = "";
			}
		}
		return vect;
	}

	double* getIdentity(int size) {
		double* matrix = new double[(int)(size)];
		for (int i = 0; i < size * size; i++) {
			matrix[i] = 0.0;
		}
		for (int i = 0; i < size; i++) {
			matrix[size * i + i] = 1.0;
		}
		return matrix;
	}
};

class funcion2 {
public:
	funcion2(double* mat_, int m_, int n_, int nnz_) {
		mat = mat_;
		m = m_;
		n = n_;
		lda = n;
		nnz = nnz_;
		info = 0;
		ja = new MKL_INT[nnz]{0};
		ia = new MKL_INT[m+1]{0};
		valuesOutCSR = new double[nnz]{0};		
		matrixOver = utils.clone(mat, m * n);
	}

	void initJOB(int * values){
		job = utils.cloneInt(values, 6);
	}

	void CSR() {		
		#pragma warning(suppress : 4996)		
		mkl_ddnscsr(job, &m, &n, matrixOver, &lda, valuesOutCSR, ja, ia, &info);
		printf("\n\tEl vector de valores es:\n");
		utils.printVectorDouble(valuesOutCSR, 1, nnz);
		printf("\n\tEl vector de los indices de columnas es:\n");
		utils.printVectorInt(ja, 1, nnz);
		printf("\n\tEl vector de suma de los elementos no-nulos(offset) es:\n");
		utils.printVectorInt(ia, 1, m+1);
	}
	
	void COO() {	
		MKL_INT * rows = new MKL_INT[nnz]{0};
		MKL_INT * cols = new MKL_INT[nnz]{0};
		double * valuesOutCOO = new double[nnz] {0};
		#pragma warning(suppress : 4996)
		mkl_dcsrcoo(job, &n, valuesOutCSR, ja, ia, &nnz, valuesOutCOO, rows, cols, &info);
		printf("\n\tEl vector de valores es:\n");
		utils.printVectorDouble(valuesOutCOO, 1, nnz);
		printf("\n\tEl vector de los indices de columnas es:\n");
		utils.printVectorInt(cols, 1, nnz);
		printf("\n\tEl vector de los indices de filas es:\n");
		utils.printVectorInt(rows, 1, nnz);
	}

	void SKY() {
		double* asky = new double[m * n]{0};
		MKL_INT* pointers = new MKL_INT[m + 1]{0};
		#pragma warning(suppress : 4996)
		mkl_dcsrsky(job, &n, valuesOutCSR, ja, ia, asky, pointers, &info);
		printf("\n\tEl vector de valores es:\n");
		utils.printVectorDouble(asky, 1, m*n);
		printf("\n\tEl vector de punteros es:\n");
		utils.printVectorInt(pointers, 1, m+1);
	}

	void CSC() {
		double* valuesOutCSC = new double[nnz] {0};
		MKL_INT* ja1 = new MKL_INT[nnz]{0};
		MKL_INT* ia1 = new MKL_INT[m+1]{0};
		#pragma warning(suppress : 4996)		
		mkl_dcsrcsc(job, &n, valuesOutCSR, ja, ia, valuesOutCSC, ja1, ia1, &info);
		printf("\n\tEl vector de valores es:\n");
		utils.printVectorDouble(valuesOutCSC, 1, nnz);
		printf("\n\tEl vector de los indices de filas es:\n");
		utils.printVectorInt(ja1, 1, nnz);
		printf("\n\tEl vector de suma de los elementos no-nulos(offset) es:\n");
		utils.printVectorInt(ia1, 1, m + 1);
	}

private:
	Utils utils;
	MKL_INT lda, m, n, info, * job, * ja, * ia, nnz;
	double * mat, * matrixOver, * valuesOutCSR;
};

int main(int argc, char* argv[]) {
	Utils utils;
	int m = 4; int n = 4;
	string matSym = 
		"1 7 0 0 "
		"0 2 8 0 "
		"5 0 3 9 "
		"0 6 0 4 ";

	double* A = utils.processMatrix(matSym, m*n);
	printf("Matriz a operar:\n");
	utils.printMatrix(A, m, n);

	int nnz = utils.countNNZ(A, m*n);
	funcion2 funcion2(A, m, n, nnz);

	printf("\nEjercicio 2_A (Resultados identicos a MATLAB):\n");

	/*
	job(0): Conversion type.
		If job(0)=0, the rectangular matrix A is converted to the CSR format;
		if job(0)=1, the rectangular matrix A is restored from the CSR format.
	job(1): index base for the rectangular matrix A.
		If job(1)=0, zero-based indexing for the rectangular matrix A is used;
		if job(1)=1, one-based indexing for the rectangular matrix A is used.
	job(2): Index base for the matrix in CSR format.
		If job(2)=0, zero-based indexing for the matrix in CSR format is used;
		if job(2)=1, one-based indexing for the matrix in CSR format is used.
	job(3): Portion of matrix.
		If job(3)=0, adns is a lower triangular part of matrix A;
		If job(3)=1, adns is an upper triangular part of matrix A;
		If job(3)=2, adns is a whole matrix A.
	job(4)=nzmax: maximum number of the non-zero elements allowed if job(0)=0.
	job(5): job indicator for conversion to CSR format.
		If job(5)=0, only array ia is generated for the output storage.
		If job(5)>0, arrays acsr, ia, ja are generated for the output storage.
	*/
	funcion2.initJOB(new int[6] {0, 0, 0, 2, nnz, 1});	
	printf("CSR");
	funcion2.CSR();

	/*
	job[0]
		If job[0]=0, the matrix in the CSR format is converted to the coordinate format;
		If job[0]=1, the matrix in the coordinate format is converted to the CSR format.
		If job[0]=2, the matrix in the coordinate format is converted to the CSR format, 
			and the column indices in CSR representation are sorted in the increasing order 
			within each row.
	job[1]
		If job[1]=0, zero-based indexing for the matrix in CSR format is used;
		If job[1]=1, one-based indexing for the matrix in CSR format is used.
	job[2]
		If job[2]=0, zero-based indexing for the matrix in coordinate format is used;
		If job[2]=1, one-based indexing for the matrix in coordinate format is used.
	job[3] -----> Resource not found
	job[4]
		job[4]=nzmax - maximum number of the non-zero elements allowed if job[0]=0.
	job[5] - job indicator.
		For conversion to the coordinate format:
			If job[5]=1, only array rowind is filled in for the output storage.
			If job[5]=2, arrays rowind, colind are filled in for the output storage.
			If job[5]=3, all arrays rowind, colind, acoo are filled in for the output storage.
		For conversion to the CSR format:
			If job[5]=0, all arrays acsr, ja, ia are filled in for the output storage.
			If job[5]=1, only array ia is filled in for the output storage.
			If job[5]=2, then it is assumed that the routine already has been called with the 
				job[5]=1, and the user allocated the required space for storing the output arrays 
				acsr and ja.
	*/
	funcion2.initJOB(new int[6]{0, 0, 0, 2, nnz, 3});
	printf("COO");
	funcion2.COO();

	/*
	job(0)
		If job(0)=0, the matrix in the CSR format is converted to the skyline format.
		If job(0)=1, the matrix in the skyline format is converted to the CSR format.
	job(1)
		If job(1)=0, zero-based indexing for the matrix in CSR format is used.
		If job(1)=1, one-based indexing for the matrix in CSR format is used.
	job(2)
		If job(2)=0, zero-based indexing for the matrix in the skyline format is used.
		If job(2)=1, one-based indexing for the matrix in the skyline format is used.
	job(3)
		For conversion to the skyline format:
			If job(3)=0, the upper part of the matrix A in the CSR format is converted.
			If job(3)=1, the lower part of the matrix A in the CSR format is converted.
		For conversion to the CSR format:
			If job(3)=0, the matrix is converted to the upper part of the matrix A in the CSR format.
			If job(3)=1, the matrix is converted to the lower part of the matrix A in the CSR format.
	job(4)
		job(4)=nzmax - maximum number of the non-zero elements of the matrix A if job(0)=0.
	job(5) - job indicator.
		Only for conversion to the skyline format:
			If job(5)=1, only arrays pointers is filled in for the output storage.
			If job(5)=0, all output arrays asky and pointers are filled in for the output storage.
	*/
	
	printf("SKY LINE");
	printf("\nMatriz inferior convertida");	
	funcion2.initJOB(new int[6]{ 0, 0, 0, 0, nnz, 0 }); //Upper matrix.
	funcion2.SKY();
	printf("\nMatriz superior convertida");
	funcion2.initJOB(new int[6]{ 0, 0, 0, 1, nnz, 0 }); //Lower matrix.
	funcion2.SKY();

	/*
	job[0]
		If job(0)=0, the matrix in the CSR format is converted to the CSC format.
		If job(0)=1, the matrix in the CSC format is converted to the CSR format.
	job[1]
		If job[1]=0, zero-based indexing for the matrix in CSR format is used;
		If job[1]=1, one-based indexing for the matrix in CSR format is used.
	job[2]
		If job[2]=0, zero-based indexing for the matrix in the CSC format is used;
		If job[2]=1, one-based indexing for the matrix in the CSC format is used.
	job[3] -----> Resource not found
	job[4] -----> Resource not found		
	job[5] - job indicator.
		For conversion to the CSC format:
			If job(5)=0, only arrays ja1, ia1 are filled in for the output storage.
			If job(5)≠0, all output arrays acsc, ja1, and ia1 are filled in for the
			output storage.
		For conversion to the CSR format:
			If job(5)=0, only arrays ja, ia are filled in for the output storage.
			If job(5)≠0, all output arrays acsr, ja, and ia are filled in for the 
			output storage.
	*/
	funcion2.initJOB(new int[6]{ 0, 0, 0, 0, 0, 1 });
	printf("CSC");
	funcion2.CSC();
	
	char a = getchar();	
	return 0;
}