#pragma once
#include "PCA.h"



PCA::PCA() {};


bool PCA::LoadData(const char* const file_name) {
	if (st != ready_to_load && st != finished) {
		AddError("LeadData error: unproper state.");
		return false;
	}
	try {
		X.ReadFromFile(file_name, TXT);
	}
	catch (std::exception& err) {
		AddError(err.what());
		return false;
	}
	st = ready_to_prepare;
	return true;
}

bool PCA::Prepare() {
	if (st == ready_to_prepare || st == finished) {
		E.Resize(X.Rows(), X.Cols());
		for (size_t j = 0; j < X.Cols(); j++) {
			double mj = Mean(X.GetCol(j));
			double sj = Dispersion(X.GetCol(j));
			for (size_t i = 0; i < X.Rows(); i++)
				X.At(i, j) = (X.At(i, j) - mj) / sj;

		}
		E = X;
		X0.Resize(X.Rows(), X.Cols());
		X0 = X;
		st = ready_for_calc;
		return true;
	}
	else {
		AddError("Prepare error: unproper state.");
		return false;
	}
}

bool PCA::Run(size_t pc) {
	if (st == ready_for_calc || st == finished) {
		for (size_t h = 0; h < pc; h++) {
			auto t = E.GetCol(h);
			ColVector<double> d(E.Cols());
			ColVector<double> p(E.Cols());

			do {

				p = Transpose(Transpose(t) * E / ScalarProduct(t, t));
				p = p / EuclidNorm(p);
				auto t_old = t;

				t = E * p / ScalarProduct(p, p);
				d = t_old - t;


			} while (EuclidNorm(d) > EPS);
			//std::cout << t * Transpose(p);
			E = E - t * Transpose(p);
			//std::cout << X;
			P = (!h ? ToMatrix(p) : Assign(P, p));
			T = (!h ? ToMatrix(t) : Assign(T, t));
		}
		st = finished;
		return true;
	}
	else {
		AddError("Prepare error: unproper state.");
		return false;
	}
}

RowVector<double> PCA::Scope() {
	if (st == finished) {
		RowVector<double> h(T.Cols());
		Matrix<double> M;
		try {
			M = Inverse(T * Transpose(T));
		}
		catch (const std::exception& e) {
			AddError(e.what());
			return RowVector<double>{};
		}
		for (size_t i = 0; i < h.Size(); i++) {
			ColVector<double> t = T.GetCol(i);
			h[i] = Transpose(t) * M * t;
		}
		return h;
	}
	else {
		AddError("Scope error: unproper state.");
		return RowVector<double>{};
	}
}

double PCA::TRV() {
	if (st == finished) {
		ColVector<double> v(E.Rows());

		for (size_t i = 0; i < v.Size(); i++)
			for (size_t j = 0; j < E.Cols(); j++)
				v[i] += E.At(i, j) * E.At(i, j);
		double v0 = Mean(v);
		return v0 / E.Cols();
	}
	else {
		AddError("TRV error: unproper state.");
		return 0;
	}
}

double PCA::ERV() {
	if (st == finished) {
		double e0 = 0.0, x0 = 0.0;
		for (size_t i = 0; i < E.Rows(); i++)
			for (size_t j = 0; j < E.Cols(); j++) {
				e0 += E.At(i, j) * E.At(i, j);
				x0 += X.At(i, j) * X.At(i, j);
			}
		return 1 - e0 / x0;
	}
	else {
		AddError("TRV error: unproper state.");
		return 0;
	}
}

RowVector<double> PCA::DispersionAV()
{
	Matrix<double> X1(T.Rows(), P.Rows()), E(X.Rows(), X.Cols());//While transposing: get_Mrows->get_Ncols
	X1 = T * Transpose(P);
	E = X0 - X1;
	RowVector<double> v(E.Rows());
	for (size_t i = 0; i < E.Rows(); i++)
	{
		double k = 0;
		for (size_t j = 0; j < E.Cols(); j++)
			k += (E.At(i,j)) * (E.At(i,j));
		v[i] = k;
	}
	return v;
}


Matrix<double> PCA::GenPcData(size_t min, size_t max) {
	if (st == ready_for_calc || st == finished) {
		Matrix<double> TE_RV(2, max - min + 1);
		for (size_t i = min; i <= max; i++) {
			Run(i);
			TE_RV.At(0, i - min) = TRV();
			TE_RV.At(1, i - min) = ERV();
		}
		return TE_RV;
	}
	else {
		AddError("GenPcData error: unproper state.");
		return Matrix<double>{};
	}
}


void PCA::FetchAllErrors(void(*call_back)(const std::string&)) {
	while (!errors.empty()) {
		call_back(errors.front());
		errors.pop();
	}
	n_err = 0;
	st = ready_to_load;
}

std::string PCA::FetchError() {
	std::string err = errors.front();
	errors.pop();
	n_err--;
	if (n_err == 0)
		st = ready_to_load;
	return err;
}

size_t PCA::ErrorsCount() {
	return n_err;
}


const Matrix<double>& PCA::GetMatrix(MatrixName m_n)const {
	switch (m_n)
	{
	case PCA::mX:
		return X;
		break;
	case PCA::mP:
		return P;
		break;
	case PCA::mT:
		return T;
		break;
	case PCA::mE:
		return E;
		break;
	default:
		break;
	}
}


void PCA::AddError(const std::string& error) {
	st == unfetched_error;
	n_err++;
	errors.push(error);
}