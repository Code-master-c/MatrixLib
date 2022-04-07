#pragma once
#include <exception>

class InvalidSize : public std::exception {
	using std::exception::exception;
};

class InvalidFileExtension: public std::exception {
	using std::exception::exception;
};

class InvalidFileStructure : public std::exception {
	using std::exception::exception;
};

class SystemExceptionWithFile: public std::exception {
	using std::exception::exception;
};