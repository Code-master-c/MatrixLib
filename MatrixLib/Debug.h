#pragma once

#include <iostream>
#include <sstream>
#include <exception>
#include <iomanip>
#include "Matrix.h"


template<typename T, typename IndexingPolicy>
std::ostream& operator<<(std::ostream& out, const Matrix<T, IndexingPolicy>& m) {
	if constexpr (std::is_same_v<IndexingPolicy, MatrixIndexingPolicy>) {
		out << "Matrix[" << m.Rows() << "x" << m.Cols() << "]";
	}
	else if constexpr (std::is_same_v<IndexingPolicy, RowVectorIndexingPolicy>) {
		out << "RowVector[" << m.Size() << "]";
	}

	else if constexpr (std::is_same_v<IndexingPolicy, ColVectorIndexingPolicy>) {
		out << "ColVector[" << m.Size() << "]";
	}
	out << " of " << typeid(T).name() << " {\n";
	for (size_t i = 0; i < m.Rows(); i++) {
		out << "\t";
		for (size_t j = 0; j < m.Cols(); j++) {
			out << std::setw(2) <<  m.At(i, j) << " ";
		}
		out << "\n";
	}

	out << "}\n";
	return out;
}