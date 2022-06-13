#include "Debug.h"
#include "PCA.h"

#define AutoValidation(object, method) \
if (!(object.method))\
object.FetchAllErrors([](const std::string& error) {std::cout << error << std::endl; });



int main() {

	PCA pca;
	AutoValidation(pca, LoadData("materials\\PCA\\data.txt"));


	/*const Matrix<double>& X = pca.GetMatrix(PCA::mX);
	const Matrix<double>& T = pca.GetMatrix(PCA::mT);
	const Matrix<double>& P = pca.GetMatrix(PCA::mP);
	const Matrix<double>& E = pca.GetMatrix(PCA::mE);*/
	//std::cout << X;

	AutoValidation(pca, Prepare());
	
	auto TR_PC = pca.GenPcData(1, 12);
	

	//std::cout << "Scope = " << pca.Scope();
	std::cout << "TRV = " << TR_PC.GetRow(0) << std::endl;
	std::cout << "PRV = " << TR_PC.GetRow(1) << std::endl;

	while (pca.ErrorsCount()) {
		std::cout << pca.FetchError() << std::endl;;
	}

	return 0;
}