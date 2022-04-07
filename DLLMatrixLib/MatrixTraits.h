#pragma once
#include "Matrix.h"

template<typename Matrix_>
struct MatrixDecay { 
	constexpr static bool isMatrix = false;
	typedef void  T;
	typedef void Policy;

	constexpr static bool isMatrixPolicy = false;
	constexpr static bool isRowVectorPolicy = false;
	constexpr static bool isColVectorPolicy = false;
	constexpr static bool isVectorPolicy = false;
};


template<typename T_, typename Policy_>
struct MatrixDecay<Matrix<T_, Policy_>> {
	constexpr static bool isMatrix = true;
	typedef  T_  T;
	typedef Policy_ Policy;

	constexpr static bool isMatrixPolicy = std::is_same_v<Policy_, MatrixIndexingPolicy>;
	constexpr static bool isRowVectorPolicy = std::is_same_v<Policy_, RowVectorIndexingPolicy>;
	constexpr static bool isColVectorPolicy = std::is_same_v<Policy_, ColVectorIndexingPolicy>;
	constexpr static bool isVectorPolicy = std::is_same_v<Policy_, RowVectorIndexingPolicy> || std::is_same_v<Policy_, ColVectorIndexingPolicy>;
};

template<typename T_>
struct MatrixDecay<Identity<T_>> {
	constexpr static bool isMatrix = true;
	using T = T_;
	using Policy =  MatrixIndexingPolicy;

	constexpr static bool isMatrixPolicy = true;
	constexpr static bool isRowVectorPolicy = false;
	constexpr static bool isColVectorPolicy = false;
	constexpr static bool isVectorPolicy = false;
};

template<typename T_>
struct MatrixDecay<Diagonal<T_>> {
	constexpr static bool isMatrix = true;
	using T = T_;
	using Policy = MatrixIndexingPolicy;

	constexpr static bool isMatrixPolicy = true;
	constexpr static bool isRowVectorPolicy = false;
	constexpr static bool isColVectorPolicy = false;
	constexpr static bool isVectorPolicy = false;
};

template<typename T_>
struct MatrixDecay<UpperTriangular<T_>> {
	constexpr static bool isMatrix = true;
	using T = T_;
	using Policy = MatrixIndexingPolicy;

	constexpr static bool isMatrixPolicy = true;
	constexpr static bool isRowVectorPolicy = false;
	constexpr static bool isColVectorPolicy = false;
	constexpr static bool isVectorPolicy = false;
};

template<typename T_>
struct MatrixDecay<LowerTriangular<T_>> {
	constexpr static bool isMatrix = true;
	using T = T_;
	using Policy = MatrixIndexingPolicy;

	constexpr static bool isMatrixPolicy = true;
	constexpr static bool isRowVectorPolicy = false;
	constexpr static bool isColVectorPolicy = false;
	constexpr static bool isVectorPolicy = false;
};

template<typename T_>
struct MatrixDecay<Symmetric<T_>> {
	constexpr static bool isMatrix = true;
	using T = T_;
	using Policy = MatrixIndexingPolicy;

	constexpr static bool isMatrixPolicy = true;
	constexpr static bool isRowVectorPolicy = false;
	constexpr static bool isColVectorPolicy = false;
	constexpr static bool isVectorPolicy = false;
};



template<typename T_>
struct RowType {
	using T = T_;
};

template<typename T_>
struct RowType<T_*> {
	using T = T_;
};

template<typename T_>
struct RowType<T_&> {
	using T = T_;
};

template<typename T_>
struct RowType<T_&&> {
	using T = T_;
};

template<typename U, typename M>
struct ChangeT {
	using T = void;
};

template<typename U, typename T_, typename P>
struct ChangeT<U, Matrix<T_, P>> {
	using T = Matrix<U, P>;
};


template<typename P, typename M>
struct ChangePolicy {
	using T = void;
};

template<typename P, typename T_, typename P_>
struct ChangePolicy<P, Matrix<T_, P_>> {
	using T = Matrix<T_, P>;
};