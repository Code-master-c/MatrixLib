#pragma once
#include <vector>
#include <type_traits>
#include "MatrixExeption.h"
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include "Policy.h"

#include "MatrixLib.h"


template<typename T, typename MatrixPolicy>
class Matrix;

template<typename T>
using RowVector = Matrix<T, RowVectorIndexingPolicy>;

template<typename T>
using ColVector = Matrix<T, ColVectorIndexingPolicy>;


template<typename T,typename IndexingPolicy = MatrixIndexingPolicy>
class Matrix
{

	static_assert(std::is_arithmetic_v<T>, "Type-argument T sould be an arithmetic type.");
	static_assert(isStandartPolicy<IndexingPolicy>::val, "Indexing policy sould be one of standart policies: MatrixIndexingPolicy, RowVectorIndexingPolicy, ColVectorIndexingPolicy");

protected:
	std::vector<T> m;
	size_t N;
	size_t M;
public:
	template<typename T, typename Policy>
	friend class Matrix;


	size_t Rows() const {
		return N;
	}

	size_t Cols() const {
		return M;
	}

	//================================================================
	// Constructors:
	// 1. Matricies:
	//  1.1 zero matrix: A() -> A[1x1],A(3) -> A[3x1], A(3, 3) -> A[3x3]
	//  1.2 copy-constructor: A<T>(B<U>), provided that type U is reduced to T
	//  1.3 init-lists: A({{...}, {...}, ...})
	//  1.4 pointer: A<T>(U* ptr, n), provided that type U is reduced to T (unsafe)
	//  1.5 vector A<T>(vector<vector<U>>), provided that type U is reduced to T (unsafe)
	// 2. Vectors
	//  2.1 zero vector: V(3) -> V[3x1] or V[1x3]
	//  2.2 copy-constructor: V<T>(V0<U>), provided that type U is reduced to T
	//  2.3 init-list: V({...})
	//  2.4 pointer: V<T>(U* ptr, n), provided that type U is reduced to T(unsafe)
	//  2.5 vector: V<T>(vector<U>), provided that type U is reduced to T (unsafe)
	//================================================================

	template<
		typename = std::enable_if_t<isMatrixPolicy<IndexingPolicy>::val, int>
	> //только для матриц.
	Matrix(size_t rows, size_t cols = 1) : N{ rows }, M{cols} {
		if (rows == 0 || cols == 0) throw InvalidSize(__FUNCSIG__" erorr: A matrix must contain at least one row and a column.");
		m.assign(rows * cols, 0);
	}

	Matrix() : N{ 1 }, M{ 1 } {
		m.assign(1, 0);
	}

	template<
		typename = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val, int>
	> //только для векторов.
	Matrix(size_t size) {
		if (size == 0)  throw InvalidSize(__FUNCSIG__" erorr: A vector must contain at least one element.");
		if constexpr (std::is_same_v<IndexingPolicy, RowVectorIndexingPolicy>) {
			N = 1; M = size;
		}
		else {
			N = size; M = 1;
		}
		m.assign(size, 0);
	}


	Matrix(const Matrix<T, IndexingPolicy>& m_) {
		N = m_.N;
		M = m_.M;
		m = m_.m;
	}

	template<typename U>
	Matrix(const Matrix<U, IndexingPolicy>& m_) {
		N = m_.N;
		M = m_.M;
		std::transform(m_.m.begin(), m_.m.end(), std::back_inserter(m), [](U a)->T {return T(a); });
	}



	template<typename U, 
		typename = std::enable_if_t<isMatrixPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int> 
	> //только для матриц.
	Matrix(const U* m_, size_t N, size_t M) : N{ N }, M{ M } {
		if (N == 0 || M == 0) throw InvalidSize(__FUNCSIG__" error: The matrix must contain at least one row and a column.");
		m.assign(m_, m_ + N * M);
	}

	template<typename U,
		typename = std::enable_if_t<isMatrixPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int>
	> //только для матриц.
	Matrix(const std::initializer_list<std::initializer_list<U>>& list) {
		N = list.size();
		if (N == 0) 
			throw InvalidSize(__FUNCSIG__" error: A matrix must contain at least one row");
		M = (*list.begin()).size();
		if (M == 0) 
			throw InvalidSize(__FUNCSIG__" error: A matrix must contain at least one column.");
		
		for (auto& it : list) {
			if (it.size() != M)
				throw std::logic_error(__FUNCSIG__" error: The lengths of the rows must match.");
			if (it.size() == 0)
				throw InvalidSize(__FUNCSIG__" error: A matrix must contain at least one row.");
		}
		m.reserve(N * M);
		for (auto& ivec : list)
			for (auto& it : ivec)
				m.push_back((T)it);
	}

	template<typename U, 
		typename  = std::enable_if_t<isMatrixPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int> 
	> //только для матриц.
	Matrix(const std::vector<std::vector<U>>& vec) {
		N = vec.size();
		if (N == 0) throw InvalidSize(__FUNCSIG__" error: A matrix must contain at least one row");
		M = vec[0].size();
		for (auto& i : vec) {
			if (i.size() != M)
				throw std::logic_error(__FUNCSIG__" error: The lengths of the strings must match.");
			if (i.size() == 0)
				throw InvalidSize(__FUNCSIG__" error: A matrix must contain at least one row.");
		}

		m.reserve(N * M);

		for (size_t i = 0; i < N; i++) for (size_t j = 0; j < M; j++)
			m.push_back(T(vec[i][j]));

	}


	template<typename U,
		typename = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val&& std::is_arithmetic_v<U>, int>
	> //только для векторов.
	Matrix(const std::initializer_list<U>& list) {
		if (list.size() == 0) 
			throw InvalidSize(__FUNCSIG__"error: A vector must contain at least one entry.");
		if constexpr (std::is_same_v<IndexingPolicy, RowVectorIndexingPolicy>) {
			N = 1; M = list.size();
		}
		else {
			N = list.size(); M = 1;
		}
		m.assign(list.begin(), list.end());
	}

	template<typename U, 
		typename  = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int>
	> //только для векторов.
	Matrix(const U* m_, size_t size) {
		if (size == 0)
			throw InvalidSize(__FUNCSIG__"error: A vector must contain at least one entry.");
		if constexpr (std::is_same_v<IndexingPolicy, RowVectorIndexingPolicy>) {
			N = 1; M = size;
		}
		else {
			N = size; M = 1;
		}
		m.assign(m_, m_ + size);
	}

	template<typename U,
		typename = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int>
	> //только для векторов.
	Matrix(const std::vector<U>& vec){
		if (vec.size() == 0)
			throw std::exception(__FUNCSIG__" erorr: The vector must contain at least one entry.");
		if constexpr (std::is_same_v<IndexingPolicy, RowVectorIndexingPolicy>) {
			N = 1; M = vec.size();
		}
		else {
			M = 1; N = vec.size();
		}
		std::transform(vec.begin(), vec.end(), std::back_inserter(m), [](U a)->T {return T(a); });
	}


	//================================================================
	// operator=
	// 1. Matrices:
	//  1.1 init-lists: A({{...}, {...}, ...})
	//		memory is not reallocated, only the maximum possible submatrix is copied.
	//  1.2 matrix: A<T>(B<U>), provided that type U is reduced to T
	//		memory is not reallocated, only the maximum possible submatrix is copied.
	//
	// 2. Vectors
	//  2.1 init-list V({...})
	//		memory is not reallocated, only the maximum possible submatrix is copied.
	//  2.2 matrix: V<T>(M<U>), provided that type U is reduced to T
	//		memory is not reallocated, only the maximum possible submatrix is copied.
	// 
	//================================================================

	template<typename U,
		typename = std::enable_if_t<isMatrixPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int>
	> //только для матриц.
	void operator=(const std::initializer_list<std::initializer_list<U>>& list) {	
		size_t i = 0, j = 0;
		for (auto& m_ : list) {
			if (i >= N) break;
			for (auto val : m_) {
				if (j >= M) break;
				m[i * N + j] = val;
				j++;
			}
			j = 0;
			i++;
		}
	}

	template<typename U,
		typename = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val && std::is_arithmetic_v<U>, int>
	> //только для векторов.
	void operator=(const std::initializer_list<U>& list) {
		size_t i = 0;
		for (auto& val : list) {
			if (i >= Size()) break;
			m[i] = val;
			i++;
		}
	}

	template<typename U, typename UIndexingPolicy>
	void operator=(const Matrix<U, UIndexingPolicy>& m_) {
		for (size_t i = 0; i < min(N, m_.N); i++) for (size_t j = 0; j < min(M, m_.M); j++)
			m[i * M + j] = T(m_.At(i, j));
	}

	void operator=(const Matrix<T, IndexingPolicy>& m_) {
		for (size_t i = 0; i < min(N, m_.N); i++) for (size_t j = 0; j < min(M, m_.M); j++)
			m[i * M + j] = m_.At(i, j);
	}

	void operator=(const Matrix<T, IndexingPolicy>&& m_) {
		N = m_.N;
		M = m_.M;
		m = m_.m;
	}


	//================================================================
	// Insdexing:
	// 1. Matricies:
	//  1.1 M<int>[1] -> RowVector<int> (!copied row!)
	// 2. Vectors:
	//  2.1: V<T>[1] -> T(&) value
	// 3. Matricies & Vectors:
	//  3.1: A<T>.At(1,2) -> A[1][2] -> T(&) value
	//================================================================

	template<typename = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val, int>> 
	T& operator[](size_t i) {
		if (i >= Size()) throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		return m[i];
	}

	template<typename =  std::enable_if_t<isVectorPolicy<IndexingPolicy>::val, int>>
	T operator[](size_t i)const {
		if (i >= Size()) throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		return m[i];
	}

	
	T& At(size_t row, size_t col) {
		if (row >= N || col >= M) throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		return m[row * M + col];
	}

	T At(size_t row, size_t col)const {
		if (row >= N || col >= M) throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		return m[row * M + col];
	}


	template<typename =  std::enable_if_t<isMatrixPolicy<IndexingPolicy>::val, int>>
	RowVector<T> operator[](size_t i)const {
		if (i >= N) throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		return RowVector<T>(std::vector<T>(m.begin() + i * M, m.begin() + i * M + M));
	}
	
	//================================================================
	// Utils
	//================================================================


	template<typename = std::enable_if_t<isVectorPolicy<IndexingPolicy>::val, int> > //только для векторов.
	size_t Size()const {
		return m.size();
	}

	void Resize(size_t rows, size_t cols) {
		if (rows == 0 || cols == 0) throw InvalidSize(__FUNCSIG__" error: The matrix must contain at least one row and a column.");
		auto m_ = *this;
		N = rows;
		M = cols;
		m.clear();
		m.resize(rows * cols);
		for (size_t i = 0; i < m_.Rows(); i++) for (size_t j = 0; j < m_.Cols(); j++)
			m[i * rows + j] = m_.At(i, j);
	}

	auto GetRow(size_t row)const {
		if (row >= N)  throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		return RowVector<T>(std::vector<T>(m.begin() + row * M, m.begin() + row * M + M));
	}

	auto GetCol(size_t col)const {
		if (col >= M)  throw std::out_of_range(__FUNCSIG__" error: index out of range.");
		ColVector<T> ret(N);
		for (size_t i = 0; i < N; i++) {
			ret[i] = m[i * M + col];
		}
		return ret;
	}

#define TXT 0
#define BINM 1
#define AUTO 2

	void ReadFromFile(const char* file_path, int format = AUTO) {
		if (format == AUTO) {
			std::string fp(file_path);

			auto i = fp.rbegin();
			while (i <= fp.rend() && *i != '.') i++; //находим позицию точки
			if (i == fp.rend() || i == fp.rbegin()) throw InvalidFileExtension(__FUNCSIG__" error: failed to recognize file extension.");

			auto ext = fp.substr(fp.rend() - (i + 1), fp.size() - 1);

			if (ext == ".txt") format = TXT;
			else if (ext == ".binm") format = BINM;
			else throw InvalidFileExtension(__FUNCSIG__" error: unknown file extension.");

		}
		if (format == TXT) {
			std::ifstream file(file_path);
			if (!file.is_open()) throw SystemExceptionWithFile(__FUNCSIG__" error: could not open the file.");
			ReadFromFStream(file, TXT);
		}
		else if (format == BINM) {
			std::ifstream file(file_path, std::ifstream::out | std::ifstream::binary);
			if (!file.is_open()) throw SystemExceptionWithFile(__FUNCSIG__" error: could not open the file.");
			ReadFromFStream(file, BINM);
		}
		else
			throw InvalidFileExtension(__FUNCSIG__" error: failed to recognize file extension.");

	}


	void ReadFromFStream(std::istream& in, int format = TXT) {
		if (format == TXT) {
			m.clear();
			N = 0;
			M = 0;
			T entry;
			std::stringstream ss;
			std::string line;

			size_t K = 0;
			while (std::getline(in, line, '\n')) {
				ss.str(line);
				while (ss >> entry) {
					m.push_back(entry);
					if (N == 0) M++;
					else K++;
				}
				if (K != M && N) throw InvalidFileStructure(__FUNCSIG__" error: rows that were not equal in length were found.");
				K = 0;
				N++;
				if constexpr (std::is_same_v<IndexingPolicy, RowVectorIndexingPolicy>) {
					if (N > 1)  throw InvalidFileStructure(__FUNCSIG__" error: the data in the file cannot be interpreted as a row vector.");
				}
				else if constexpr (std::is_same_v<IndexingPolicy, ColVectorIndexingPolicy>) {
					if (M > 1)  throw InvalidFileStructure(__FUNCSIG__" error: the data in the file cannot be interpreted as a col vector.");
				}
				ss.clear();
			}
			if (N == 0 || M == 0) throw InvalidFileStructure(__FUNCSIG__" error: the data in the file cannot be interpreted as a matrix object.");
		}
		else if (format == BINM) {
			in.read((char*)&N, sizeof(size_t));
			in.read((char*)&M, sizeof(size_t));
			if (N == 0 || M == 0) throw InvalidFileStructure(__FUNCSIG__" error: the data in the file cannot be interpreted as a matrix object.");
			m.clear();
			m.reserve(M * N);
			double entry;

			size_t k = 0;

			while (in.read((char*)&entry, sizeof(double)) && k < N * M) {
				m.push_back(entry);
				k++;
			}
			if (k != N * M) throw InvalidFileStructure(__FUNCSIG__" error: there were fewer elements predicted by the file description than there actually are. The file may have been corrupted.");
		}
		else
			throw std::logic_error(__FUNCSIG__" error: unknown data format.");
	}

	void Save(std::ostream& out, bool format = TXT)const{
		if (format == TXT) {
			for (size_t i = 0; i < N; i++) {
				for (size_t j = 0; j < M; j++)
					out << m[i * M + j] << '\t';
				out << '\n';
			}
			out.flush();
		}
		else if (format == BINM) {
			int description = sizeof(T);
			out.write((char*)&N, sizeof(size_t));
			out.write((char*)&M, sizeof(size_t));
			for (size_t i = 0; i < N; i++) for (size_t j = 0; j < M; j++) {
				double entry = m[i * M + j];
				out.write((char*)&entry, sizeof(double));
			}
			out.flush();
		}
		else
			throw std::exception(__FUNCSIG__" error: unknown data format.");
	}

	void Save(const char* file_path, int format = AUTO)const {
		if (format == AUTO) {
			std::string fp(file_path);
			auto i = fp.rbegin();
			while (i <= fp.rend() && *i != '.') i++; //находим позицию точки
			if (i == fp.rend() || i == fp.rbegin()) InvalidFileExtension(__FUNCSIG__" error: failed to recognize file extension.");

			auto ext = fp.substr(fp.rend() - (i + 1), fp.size() - 1);

			if (ext == ".txt") format = TXT;
			else if (ext == ".binm") format = BINM;
			else throw InvalidFileExtension(__FUNCSIG__" error: unknown file extension.");

		}

		if (format == TXT) {
			std::ofstream file(file_path);
			if (!file.is_open()) throw SystemExceptionWithFile(__FUNCSIG__" error: could not open the file.");
			Save(file, TXT);
			file.close();
		}
		else if (format == BINM) {
			std::ofstream file(file_path, std::ofstream::out | std::ofstream::binary);
			if (!file.is_open()) throw SystemExceptionWithFile(__FUNCSIG__" error: could not open the file.");
			Save(file, BINM);
			file.close();
		}
		else
			throw SystemExceptionWithFile(__FUNCSIG__" error: unknown file extension.");
	}

};


struct EntryPosition {
	size_t row;
	size_t col;
};


template<typename T>
class Identity :public Matrix<T> {
public:
	Identity(size_t n) : Matrix<T>(n,n) {
		for (size_t i = 0; i < n; i++) {
			this->m[i * this->M + i] = 1;
		}
	}
};


template<typename T>
class Diagonal :public Matrix<T> {
public:
	Diagonal(T a, size_t n) : Matrix<T>(n, n) {
		for (size_t i = 0; i < n; i++) {
			this->m[i * this->M + i] = a;
		}
	}
};


//верхняя и нижняя треугольные матрицы, симметричная матрица

template<typename T>
class UpperTriangular :public Matrix<T> {
public:
	UpperTriangular(size_t n, const std::vector<T>& entries) : Matrix<T>(n, n) {		
		if (entries.size() < ((1 + n) * n) / 2) throw InvalidSize(__FUNCSIG__" error: not enough elements were transferred to form the matrix.");
		size_t k = 0;
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < n; j++)
				if (j < i)
					this->m[i * this->N + j] = 0;
				else {
					this->m[i * this->N + j] = entries[k];
					k++;
				}
	}
};

template<typename T>
class LowerTriangular :public Matrix<T> {
public:
	LowerTriangular(size_t n, const std::vector<T>& entries) : Matrix<T>(n, n) {
		size_t k = 0;
		if (entries.size() < ((1 + n) * n) / 2) throw InvalidSize(__FUNCSIG__" error: not enough elements were transferred to form the matrix.");
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < n; j++)
				if (j > i)
					this->m[i * this->N + j] = 0;
				else {
					this->m[i * this->N + j] = entries[k];
					k++;
				}
	}
};

template<typename T>
class Symmetric :public Matrix<T> {
public:
	Symmetric(size_t n, const std::vector<T>& entries) : Matrix<T>(n, n) {
		/*
		* 3 4 1 5 4 7 ->

		   3 4 1
			 5 4
			   7
		*/
		if (entries.size() < ((1 + n) * n) / 2) throw InvalidSize(__FUNCSIG__" error: not enough elements were transferred to form the matrix.");
		size_t k = 0;
		for (size_t i = 0; i < n; i++)
			for (size_t j = i; j < n; j++) {
				this->m[i * this->N + j] = entries[k];
				k++;
			}

		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < i; j++)
				this->m[i * this->N + j] = this->m[j * this->N + i];
	}
};