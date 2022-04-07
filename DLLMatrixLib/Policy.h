#pragma once



struct MatrixIndexingPolicy {};
struct RowVectorIndexingPolicy {};
struct ColVectorIndexingPolicy {};

template<typename T>
struct isStandartPolicy {
	constexpr bool static val = false;
};

template<>
struct isStandartPolicy<MatrixIndexingPolicy> {
	constexpr bool static val = true;
};


template<>
struct isStandartPolicy<RowVectorIndexingPolicy> {
	constexpr bool static val = true;
};

template<>
struct isStandartPolicy<ColVectorIndexingPolicy> {
	constexpr bool static val = true;
};

template<typename Policy>
struct isVectorPolicy {
	constexpr bool static val = std::is_same_v<Policy, RowVectorIndexingPolicy> || std::is_same_v<Policy, ColVectorIndexingPolicy>;
};


template<typename Policy>
struct isMatrixPolicy {
	constexpr bool static val = std::is_same_v<Policy, MatrixIndexingPolicy>;
};

