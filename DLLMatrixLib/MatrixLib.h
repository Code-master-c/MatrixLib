#pragma once


#ifdef DLLMATRIXLIB_EXPORTS
#define MATRIXLIB_API __declspec(dllexport)
#else
#define MATRIXLIB_API __declspec(dllimport)
#endif

#include "Matrix.h"
#include "MatrixUtility.h"
#include "Debug.h"


extern "C" template class MATRIXLIB_API Matrix<double>;

extern "C" MATRIXLIB_API void foo();