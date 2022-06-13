#include "Debug.h"
#include "PCA.h"

#define AutoValidation(object, method) \
if (!(object.method))\
object.FetchAllErrors([](const std::string& error) {std::cout << error << std::endl; });

int main() {

	PCA pca;
	AutoValidation(pca, LoadData("materials\\PCA\\data.txt"));


	const Matrix<double>& X = pca.GetMatrix(PCA::mX);
	const Matrix<double>& T = pca.GetMatrix(PCA::mT);
	const Matrix<double>& P = pca.GetMatrix(PCA::mP);
	const Matrix<double>& E = pca.GetMatrix(PCA::mE);
	std::cout << X;
	
	AutoValidation(pca, Prepare());
	AutoValidation(pca, Run());
	std::cout << "X_scaled = " << X;
	std::cout << "X = T Pt + E = " << T * Transpose(P) + E;
	std::cout << "E = " << E;
	std::cout << "T = " << T;
	std::cout << "P = " << P;
	std::cout << "Scope = " << pca.Scope();
	std::cout << "ERV = " << pca.ERV() << std::endl;
	std::cout << "TRV = " << pca.TRV() << std::endl;
	std::cout << "DispersionAV = " << pca.DispersionAV() << std::endl;


	while (pca.ErrorsCount()) {
		std::cout << pca.FetchError() << std::endl;;
	}

	return 0;
}