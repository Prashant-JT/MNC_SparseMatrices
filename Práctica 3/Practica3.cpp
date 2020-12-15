#include <cstdio>
#include <cstdlib>
#include <string>
#include <random>
#include <ctime>
#include <fstream>

#include <mkl.h>

using namespace std;

class Utils {
public:
	void printMatrix(double* mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (mat[i] <= 0) {
				printf("%8.3f  ", mat[i]);
			}
			else {
				printf("%8.2f  ", mat[i]);
			}
			if ((i + 1) % n == 0) {
				printf("\n\t\t");
			}
		}
		printf("\n");
	}

	void printVectorInt(int* mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (i == (m * n) - 1) {
				printf("%d\n", mat[i]);
			}
			else {
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

	int countNNZ(double* mat, int size) {
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

	double* genSparseMatrix(int m, int n, int vMin, int vMax, int zDensity) {
		if (m <= 1 || n <= 1) {
			m = 2; n = 2;
		}
		if (zDensity >= m * n) zDensity = (m * n) - 1;
		if (zDensity <= 0) zDensity = 1;		
		if (vMin <= 0 || vMax <= vMin) {
			vMin = 1; vMax = 10;
		}

		std::default_random_engine generator;
		std::uniform_real_distribution<double> aleatorio(vMin, vMax);

		double* sparseMat = new double[m * n]{ 0 };
		while (zDensity > 0) {
			int pos = (rand() % (m*n));
			if (sparseMat[pos] == 0.0) {
				sparseMat[pos] = aleatorio(generator);
				zDensity--;
			}
		}
		return sparseMat;
	}

	void saveFile(int size, double* vectBlas, double* vectBlasCSR, double* vectBlasCOO, string name) {
		ofstream myfile(name);
		if (myfile.is_open())
		{
			for (int x = 0; x < size; x++) {
				myfile << vectBlas[x] << ", " << vectBlasCSR[x] << ", " << vectBlasCOO[x] << "\n";
			}
			myfile.close();
		}
	}
};

class funcion3 {
public:
	funcion3(int timesCont) {
		//Measure times
		cont = 0;
		timesBlas = new double[timesCont];
		timesSparseBlasCSR = new double[timesCont];
		timesSparseBlasCOO = new double[timesCont];
	}
	void setValues(double * matA, double* matB, int m_, int n_, int nnzA_, int nnzB_) {		
		matrixA = matA;
		matrixB = matB;
		m = m_;
		n = n_;		
		lda = n;
		nnzA = nnzA_;
		nnzB = nnzB_;
		info = 0;
		//Init Matrix A CSR
		jobA = new MKL_INT[6]{ 0, 0, 0, 2, nnzA, 1 };
		jaA = new MKL_INT[nnzA]{ 0 };
		iaA = new MKL_INT[m + 1]{ 0 };
		valuesOutCSRA = new double[nnzA] {0};
		//Init Matrix B CSR
		jobB = new MKL_INT[6]{ 0, 0, 0, 2, nnzB, 1 };
		jaB = new MKL_INT[nnzB]{ 0 };
		iaB = new MKL_INT[m + 1]{ 0 };
		valuesOutCSRB = new double[nnzB] {0};	
		//Init Matrix A COO 
		rowsA = new MKL_INT[nnzA]{ 0 };
		colsA = new MKL_INT[nnzA]{ 0 };
		jobACOO = new int[6]{ 0, 0, 0, 2, nnzA, 3 };
		valuesOutCOOA = new double[nnzA] {0};
		//Init Matrix B COO
		rowsB = new MKL_INT[nnzB]{ 0 };
		colsB = new MKL_INT[nnzB]{ 0 };
		jobBCOO = new int[6]{ 0, 0, 0, 2, nnzB, 3 };
		valuesOutCOOB = new double[nnzB] {0};
	}

	void setCont() {
		cont++;
	}

	void CSR() {
		//Fill three vectors of matrix A.
		#pragma warning(suppress : 4996)		
		mkl_ddnscsr(jobA, &m, &n, matrixA, &lda, valuesOutCSRA, jaA, iaA, &info);
		
		//Fill three vectors of matrix B.
		#pragma warning(suppress : 4996)		
		mkl_ddnscsr(jobB, &m, &n, matrixB, &lda, valuesOutCSRB, jaB, iaB, &info);
	}

	void COO() {		
		#pragma warning(suppress : 4996)
		mkl_dcsrcoo(jobACOO, &n, valuesOutCSRA, jaA, iaA, &nnzA, valuesOutCOOA, rowsA, colsA, &info);

		#pragma warning(suppress : 4996)		
		mkl_dcsrcoo(jobBCOO, &n, valuesOutCSRB, jaB, iaB, &nnzB, valuesOutCOOB, rowsB, colsB, &info);
	}
	
	void blas(bool change) {
		MKL_INT k, alpha, beta, ldb, ldc;
		CBLAS_LAYOUT layout;
		CBLAS_TRANSPOSE trans;
		layout = CblasRowMajor;
		trans = CblasNoTrans;
		k = m;
		alpha = 1;
		beta = 0;
		lda = k;
		ldb = n;
		ldc = n;
		double* matrixC = new double[m*n]{ 0 };
		
		double start, fin = dsecnd();

		if (change) { //Opera A*B
			printf("Calculo de A*B usando BLAS (sin compresion):\n");
			start = dsecnd();
			cblas_dgemm(layout, trans, trans, m, n, k, alpha, matrixA, lda, matrixB, ldb, beta, matrixC, ldc);
			fin = (dsecnd() - start);
		} else { //Opera B*A
			printf("Calculo de B*A usando BLAS (sin compresion):\n");
			start = dsecnd();
			cblas_dgemm(layout, trans, trans, m, n, k, alpha, matrixB, ldb, matrixA, lda, beta, matrixC, ldc);
			fin = (dsecnd() - start);
		}	
		timesBlas[cont] = fin;
		printf("\t--(CBLAS) El tiempo empleado para A*B (%dx%d) con NNZeros: %d fue: %f s.\n", m, n, nnzA, fin);
		//utils.printMatrix(matrixC, m, n);
	}

	void sparseBlasCSR(bool change) {
		MKL_INT k, ldb, ldc;
		k = m;
		ldb = n;
		ldc = k;
		char trans = 'N'; //No tranponse
		char* matdescra = new char[6]{'G', '_', 'N', 'C', '_', '_'};
		double alpha = 1.0;
		double beta = 0.0;
		double* matrixC = new double[m * n]{ 0.0 };
				
		double start, fin = dsecnd();

		if (change){ //Opera A*B con A comprimida y B sin comprimir.
			printf("Calculo de A*B usando Sparse BLAS (solo A comprimida):\n");
			start = dsecnd();
			#pragma warning(suppress : 4996)	
			mkl_dcsrmm(&trans, &m, &n, &k, &alpha, matdescra, valuesOutCSRA, jaA, iaA, &(iaA[1]), matrixB, &ldb, &beta, matrixC, &ldc);
			fin = (dsecnd() - start);
		} else { //Opera B*A con A sin comprimir y B comprimida.
			printf("Calculo de B*A usando Sparse BLAS (solo B comprimida):\n");
			start = dsecnd();
			#pragma warning(suppress : 4996)	
			mkl_dcsrmm(&trans, &m, &n, &k, &alpha, matdescra, valuesOutCSRB, jaB, iaB, &(iaB[1]), matrixA, &lda, &beta, matrixC, &ldc);
			fin = (dsecnd() - start);
		}		
		timesSparseBlasCSR[cont] = fin;
		printf("\t--(CSR) El tiempo empleado para A*B (%dx%d) con NNZeros: %d fue: %f s.\n", m, n, nnzA, fin);
		//utils.printMatrix(matrixC, m, n);
	}

	void sparseBlasCOO(bool change) {
		MKL_INT k, ldb, ldc;
		k = m;
		ldb = n;
		ldc = k;
		char trans = 'N'; //No tranponse
		char* matdescra = new char[6]{ 'G', '_', 'N', 'C', '_', '_' };
		double alpha = 1.0;
		double beta = 0.0;
		double* matrixC = new double[m * n]{ 0.0 };

		double start, fin = dsecnd();

		if (change) { //Opera A*B con A comprimida y B sin comprimir.
			printf("Calculo de A*B usando Sparse BLAS (solo A comprimida):\n");
			start = dsecnd();
			#pragma warning(suppress : 4996)			
			mkl_dcoomm(&trans, &m, &n, &k, &alpha, matdescra, valuesOutCOOA, rowsA, colsA, &nnzA, matrixB, &ldb, &beta, matrixC, &ldc);
			fin = (dsecnd() - start);
		}
		else { //Opera B*A con A sin comprimir y B comprimida.
			printf("Calculo de B*A usando Sparse BLAS (solo B comprimida):\n");
			start = dsecnd();
			#pragma warning(suppress : 4996)	
			mkl_dcoomm(&trans, &m, &n, &k, &alpha, matdescra, valuesOutCOOB, rowsB, colsB, &nnzA, matrixA, &lda, &beta, matrixC, &ldc);
			fin = (dsecnd() - start);
		}
		timesSparseBlasCOO[cont] = fin;
		printf("\t--(COO) El tiempo empleado para A*B (%dx%d) con NNZeros: %d fue: %f s.\n", m, n, nnzA, fin);
		//utils.printMatrix(matrixC, m, n);
	}

	double* getTimes(int blas) {
		if (blas == 0) return timesBlas;
		if (blas == 1) return timesSparseBlasCSR;
		return timesSparseBlasCOO;
	}

private:
	Utils utils;
	MKL_INT lda, m, n, info;
	MKL_INT * jobA, * jaA, * iaA, nnzA;
	MKL_INT * jobB, * jaB, * iaB, nnzB;
	MKL_INT * jobACOO, * rowsA, * colsA;
	MKL_INT * jobBCOO, * rowsB, * colsB;
	double * matrixA, * matrixB, * valuesOutCSRA, * valuesOutCSRB, * valuesOutCOOA, * valuesOutCOOB, * timesBlas, * timesSparseBlasCSR, * timesSparseBlasCOO;
	int cont;
};

int main(int argc, char* argv[]) {
	srand((unsigned)time(NULL));
	int min = 1; int max = 30000; int jump = 100;
	Utils utils;
	int sizeVects = (max / jump);
	funcion3 funcion3(sizeVects); //Inicializa vectores de tiempos.

	for (int x = (max*min); x > min; x = x-jump) {
		int m = 1000; int n = 1000;

		int zDensityA = x; int vMinA = 10; int vMaxA = 100;
		double* matA = utils.genSparseMatrix(m, n, vMinA, vMaxA, zDensityA);
		int zDensityB = x; int vMinB = 33; int vMaxB = 177;
		double* matB = utils.genSparseMatrix(m, n, vMinB, vMaxB, zDensityB);
		
		/*
		printf("Matriz A generada (%dx%d), con NNZeros=%d:\n", m, n, zDensityA);
		utils.printMatrix(matA, m, n);

		printf("Matriz B generada (%dx%d), con NNZeros=%d:\n", m, n, zDensityB);
		utils.printMatrix(matB, m, n);
		*/

		int nnzA = utils.countNNZ(matA, m * n);
		int nnzB = utils.countNNZ(matB, m * n);
		funcion3.setValues(matA, matB, m, n, nnzA, nnzB);		

		bool change = true; //Intercambia A*B a B*A (Entonces B es comprimida en lugar de A)
		funcion3.blas(change);
		funcion3.CSR(); //Comprime las matrices disperas para realizar las operaciones (CSR).
		funcion3.sparseBlasCSR(change);
		funcion3.COO(); //Comprime las matrices disperas para realizar las operaciones (COO).
		funcion3.sparseBlasCOO(change);
		printf("--------------------------------------------\n");
		funcion3.setCont();
	}

	utils.saveFile(sizeVects, funcion3.getTimes(0), funcion3.getTimes(1), funcion3.getTimes(2), "P3_times_Blas-CSR-COO.csv");

	char a = getchar();
	return 0;
}