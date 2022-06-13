#pragma once
#include "Matrix.h"
#include "MatrixUtility.h"
#include "Debug.h"
#include <queue>

#define EPS 1e-8

class PCA {
public:

	PCA();

	enum State {
		ready_to_load,
		ready_to_prepare,
		ready_for_calc,
		unfetched_error,
		finished
	};

	enum MatrixName {
		mX,mP,mT, mE
	};

	bool LoadData(const char* const file_name);

	bool Prepare();

	bool Run(size_t pc = 12);

	RowVector<double> Scope();

	double TRV();

	double ERV();

	Matrix<double> GenPcData(size_t min, size_t max);

	RowVector<double> DispersionAV();

	void FetchAllErrors(void(*call_back)(const std::string&));

	std::string FetchError();

	size_t ErrorsCount();

	const Matrix<double>& GetMatrix(MatrixName m_n)const;

private:
	void AddError(const std::string& error);

	State st = ready_to_load;
	size_t n_err;
	std::queue<std::string> errors;

	Matrix<double> X;
	Matrix<double> X0;
	Matrix<double> T;
	Matrix<double> P;
	Matrix<double> E;
};