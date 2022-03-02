/*
На основе динамических массивов, типа std::vector или другого стандартного контейнера реализовать типы данных для представления матриц и векторов произвольной размерности.

Создать класс Матрица. (Yes)
Вектор рассматривать как матрицу размера 1xN или Nx1.(Yes)

написать подклассы для матриц указанных видов - Единичная, диагональная, верхняя и нижняя треугольные матрицы, симметричная матрица. (Yes)

Использовать перегрузку операторов для реализации операций сложения, вычитания и умножения векторов и матриц. (Yes)

При оформлении библиотеки классов вынести определения классов и реализацию в отдельные файлы (h и cpp).

Использовать механизм исключений для обработки нештатных ситуаций.(Yes)

Реализовать операции над векторами и матрицами: (Yes)

Сложение матриц, вычитание матриц и умножение матриц на число (Yes)
Умножение матриц (Yes)
произведение Адамара (Yes)

перегрузку оператор вывода в поток (<<) (YEs)
*/

/*
Дополнить библиотеку классов, разработанную в рамках ЛР1. Реализовать операции над векторами и матрицами:
След матрицы (Yes)

Определитель матрицы (методом Гаусса) (Yes)

Скалярное произведение векторов (Yes)

Норма вектора (евклидова норма, максимальная норма) (Yes)

Норма матрицы (норма Фробениуса) (Yes)
*/

/*
Дополнить библиотеку классов, разработанную в рамках ЛР1 и ЛР2. Реализовать операции над векторами и матрицами:

Угол между векторами (Yes)

Ранг Матрицы (с помощью алгоритма Гаусса) (Yes)

Обратная матрица (если существует) (Yes)

Транспонирование матрицы (Yes)
*/

/*
Дополнить библиотеку классов, разработанную в рамках ЛР1-ЛР3.
Реализовать функциональность для загрузки матриц и векторов из файла и сохранении матриц и векторов в файл. (Yes)
Реализовать 2 режима работы: чтение/запись в текстовый и бинарный файл. (Yes)
Для текстового режима использовать перегрузку операторов >> и <<. (Yes)
Для бинарного режима реализовать в классах дополнительные методы. (Yes)
Использовать механизм исключений для обработки нештатных ситуаций при работе с файлами. (Yes)
*/

#include <iostream>
#include "Debug.h"
#include "Matrix.h"
#include "MatrixUtility.h"

const char* const outline = "=====================================================================";

void subcounter(const char* str, bool reset = false) {
	static char c = 'a';
	if (reset) { c = 'a'; return; }
	std::cout << '\t' << c << ')' << str << std::endl;
	c++;
}

void counter(const char* str, bool reset = false) {
	static int c = 1;
	if (reset) { c = 1; return; }
	std::cout << c << ") " << str << std::endl;
	subcounter("", true);
	c++;
}

void begin_section(const char* title) {
	std::cout << outline << std::endl;
	std::cout << '=' << std::setw(strlen(outline) - 2)<<std::left << title << '=' << std::endl;
	std::cout << outline << std::endl << std::endl;
	counter("", true);
}

void end_section() {
	std::cout << outline << std::endl << std::endl;
}


 

int main() {
	setlocale(LC_ALL, "rus");
	//=====================================================
	// тест конструкторов по умолчанию
	//=====================================================
	begin_section("Конструкторы матричных объектов по умолчанию");
	{
		counter("Конструктор матрицы по умолчанию: Matrix<int> m");
		Matrix<int> m; //матрица из int размера 1x1, заполненная нулями
		std::cout << m;

		counter("Конструктор вектора-строки по умолчанию: RowVector<float> rv");
		RowVector<float> rv; //вектор-строка размера 1, заполненный нулями
		std::cout << rv;

		counter("Конструктор вектора-столбца по умолчанию: RowVector<float> cv");
		RowVector<float> cv; //вектор-строка размера 1, заполненный нулями
		std::cout << cv;
	}
	end_section();
	//=====================================================

	//=====================================================
	// конструирование матричных объектов по их размеру.
	//=====================================================
	begin_section("Конструкторы матричных объектов с двумя аргументами");
	{
		counter("Конструктор матрицы с двумя аргументами: Matrix<int> m(3,3)");
		Matrix<int> m1(3, 3); //матрица из int размера 3x3, заполненная нулями
		std::cout << m1;

		counter("Конструктор матрицы с одним аргументом: Matrix<int> m(3) = Matrix<int> m(3,1)");
		Matrix<int> m2(3); //матрица из int размера 3x1, заполненная нулями
		std::cout << m2;

		counter("Конструктор вектора-столбца с одним аргументом: ColVector<float> cv(3)");
		ColVector<float> cv(3); //вектора-столбец, размера 3, заполненный нулями
		std::cout << cv;

		counter("Конструктор вектора-строки с одним аргументом: RowVector<float> rv(3)");
		RowVector<float> rv(3); //вектора-строки, размера 3, заполненный нулями
		std::cout << rv;
	}
	end_section();
	//=====================================================

	//=====================================================
	// конструкторы копирования
	//=====================================================
	begin_section("Конструкторы копирования матричных объектов");
	{
		counter("Конструктор копирования по умолчанию для матриц: Matrix<double> a = Matrix<double>{3,3}");
		Matrix<double> a = Matrix<int>{ 3,3 };
		std::cout << a;

		counter("Конструктор копирования матрицы с другим типом элементов: Matrix<int> b = Matrix<double>{ 2,3 }");
		Matrix<int> b = Matrix<float>{ 2,3 };
		std::cout << b;

		counter("Конструктор копирования по умолчанию для векторов: RowVector<int> rv = RowVector<float>{ 4 }");
		RowVector<int> rv = RowVector<double>{ 4 };
		std::cout << rv;

		//RowVector<int> rv = ColVector<float>{ 4 }; //ошибка!
		//std::cout << rv;

	}
	end_section();
	//=====================================================

	//=====================================================
	// конструкторы c C-массивом в качестве аргумента
	//=====================================================
	begin_section("Конструкторы матричных элементов c C-массивом в качестве аргумента");
	{
		int m_[] = { 1,2,3,4 };

		counter("Создание матрицы из элементов массива: Matrix<float> m(m_, 2, 2), где int m_[] = { 1,2,3,4 };");
		Matrix<float> m(m_, 2, 2);
		std::cout << m;

		counter("Создание вектороа из элементов массива: RowVector<int> v(m_, 3);");
		RowVector<int> v(m_, 3);
		std::cout << v;

	}
	end_section();
	//=====================================================

	//=====================================================
	// конструкторы c вектором в качестве аргумента
	//=====================================================
	begin_section("Конструкторы матричных элементов c C-массивом в качестве аргумента");
	{

		counter("Создание матрицы из элементов вектора из векторов: Matrix<float> m(m_), где std::vector<std::vector<int>> m_ = { {1,2},{2,3} };");
		std::vector<std::vector<int>> m_ = { {1,2, 0},{2,3, 1} };
		Matrix<float> m(m_);
		std::cout << m;

		subcounter("Выбрасывает исключение, если вектор не описывает матрицу: Matrix<float> m(incorrect_m_), где std::vector<std::vector<int>> incorrect_m_ = { {1,2},{2,3,1} }.");
		std::vector<std::vector<int>> incorrect_m_ = { {1,2},{2,3,1} };
		try {
			Matrix<float> m(incorrect_m_);
		}
		catch (std::exception& e) {
			std::cout << e.what() << std::endl;
		}

		counter("Создание вектора из одномерного вектора: RowVector<double> v(v_), где std::vector<float> v_ = { 1.0, 0.2, 3 }.");
		std::vector<float> v_ = { 1.0, 0.2, 3 };
		RowVector<double> v(v_);
		std::cout << v;
	}
	end_section();
	//=====================================================

	//=====================================================
	// конструкторы cо списком инициализации.
	//=====================================================
	begin_section("Конструкторы матричных элементов cо списком инициализации.");
	{
		counter("Создание матрицы от двумерного списка инициализации.");
		Matrix<float> m = {
			{1.2, 2.2, 3.2},
			{1.2, 0.3, 0.5},
			{9.8, 0.3, 0.9}
		};
		std::cout << m;

		counter("Воздание веттора от одномерного списка инициализации.");
		RowVector<int> v = { 1,2,3 };
		std::cout << v;

	}
	end_section();
	//=====================================================

	//=====================================================
	// Функции аксессоры
	//=====================================================
	begin_section("Функции аксессоры.");
	{
		counter("Функции матриц: .Rows(), .Cols()");
		Matrix<float> m = {
			{1.2, 2.2, 3.2},
			{1.2, 0.3, 0.5},
			{9.8, 0.3, 0.9}
		};
		std::cout << m << std::endl;
		std::cout << "Rows: " << m.Rows() << ", Cols: " << m.Cols() << std::endl;

		counter("Функции векторов: .Rows(), .Cols(), .Size()");
		RowVector<float> v = m[1];

		std::cout << v << std::endl;
		std::cout << "Rows: " << v.Rows() << ", Cols: " << v.Cols() << ", Size: " << v.Size() << std::endl;

		counter("Доступ к элементам с помощю функции At. Изменение элементов: m2.At(1, 1) = 0;");
		Matrix<float> m2 = {
			{1, 2, 3},
			{3, 5, 6},
			{7, 8, 9}
		};
		std::cout << m2 << std::endl;
		m2.At(1, 1) = 0;
		std::cout << m2 << std::endl;

		counter("Функция At для векторов: v2.At(0, 3) = 0.");
		RowVector<int> v2 = { 1,2,3,4 };
		std::cout << v2 << std::endl;
		v2.At(0, 3) = 0;
		std::cout << v2 << std::endl;

		counter("operator[] для матриц: m3[1]");
		Matrix<float> m3 = {
			{1, 2, 3},
			{3, 5, 6},
			{7, 8, 9}
		};
		std::cout << m3 << std::endl;
		std::cout << "Вторая строка: " << m3[1] << std::endl;


		counter("Функции .GetRow(), .GetCol(): std::cout << m3.GetRow(1) << m3.GetCol(1);");
		std::cout << m3.GetRow(1) << m3.GetCol(1);

		counter("operator[] для векторов: v3[0] = 0.");
		RowVector<int> v3 = { 1,2,3,4 };
		std::cout << v3 << std::endl;
		v3[1] = 0;
		std::cout << v3 << std::endl;

	}
	end_section();
	//=====================================================

	//=====================================================
	// Функии сохранения в файл, загрузки из файла. Сохранение и считывание из матричного потока.
	//=====================================================
	begin_section("Функии сохранения в файл, загрузки из файла. Сохранение и считывание матриц из потока.");
	{
		counter("Функции .Save и .ReadFromFile для текстовых и бинарных файлов");
		Matrix<float> m = {
			{1.2, 2.2, 3.2},
			{1.2, 0.3, 0.5},
			{9.8, 0.3, 0.9}
		};
		m.Save("materials/m_save.txt"); //текстовый формат
		m.Save("materials/m_save.binm");//бинарный формат

		Matrix<double> b;
		b.ReadFromFile("materials/m_save.txt");
		std::cout << "loaded from .txt: " << b << std::endl;
		b.ReadFromFile("materials/m_save.binm");
		std::cout << "loaded from .binm" << b << std::endl;

		counter("Сохранение и считывание из матриц из потока.");
		std::fstream out("materials/stream.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
		if (!out.is_open()) return 1;

		Matrix<float> outM =
		{ {1,2,3},
			{3,3,4},
			{2,3,1}
		};

		Matrix<int> inM;

		out << outM;
		out.seekg(0);
		out >> inM;
		std::cout << inM;
		out.close();

		counter("Вызов функций .Save() и ..ReadFromFile с дополнительным параметром.");
		m.Save("materials/m_save.strage_ext1", TXT);
		m.Save("materials/m_save.strage_ext2", BINM);

		m.ReadFromFile("materials/m_save.strage_ext1", TXT);
		std::cout << m << std::endl;

		m.ReadFromFile("materials/m_save.strage_ext2", BINM);
		std::cout << m << std::endl;

	}
	end_section();
	//=====================================================

	//=====================================================
	// Специальные виды матриц.
	//=====================================================
	begin_section("Специальные виды матриц");
	{
		counter("Единичная матрица размера n");
		std::cout << Identity<int>(10);

		counter("Верхне-треуголная матрица.");
		std::vector<int> m1_ = 
		{1, 2, 3,
		    1, 2,
		       2};

		std::cout << UpperTriangular<int>(3, m1_);

		counter("Нижне-треуголная матрица.");
		std::vector<int> m2_ =
		{ 1,
		  1, 2,
		  2, 9 ,1 };

		std::cout << LowerTriangular<int>(3, m2_);

		counter("Симетричная матрица");
		std::vector<int> m3_ =
		{
			    1,     2,   3,
			  /*2*/    2,   1,
			  /*3*/  /*1*/  0 
		};

		std::cout << Symmetric(3, m3_);
	}
	end_section();
	//=====================================================
	
	//=====================================================
	// Работа с матрицами.
	//=====================================================
	begin_section("Работа с матрицами");
	{
		counter("Умножение матриц.");
		{
			Matrix<float> a = { {1.0,2.0,3.0}, {0.2, 1.3, 0.4}, {0.0 ,0.0, 1.0} };

			std::cout << "a = " << a;

			Matrix<double> b = { {4.0,2.0}, {2.2, 1.3}, {2.0 ,0.0} };

			std::cout << "b = " << b << std::endl;

			std::cout << "a * b = " << a * b << std::endl;

			ColVector<int> v = { 1,2,3 };

			std::cout << "c = " << v << std::endl;

			std::cout << "a * v = " << a * v << std::endl;
		}

		counter("Сложение, вычитание матриц");
		{
			Matrix<int> A = { {1,2}, {3,4} };
			Identity<int> B(2);

			std::cout << "A = " << A <<std::endl;
			std::cout << "B = " << B << std::endl;

			std::cout << "A + B = " << A + B << std::endl;
			std::cout << "A - B = " << A - B << std::endl;

		}

		counter("Сложение, вычитание векторов");
		{
			ColVector<int> A = { 1,2,3 };

			std::cout << "A = " << A << std::endl;

			std::cout << "A + 0.5 * A = " << A + 0.5 * A << std::endl;
			std::cout << "A - 0.3 * A = " << A - 0.3 * A << std::endl;
		}

		counter("Create, ToMatrix, ScalarProduct, Transpose, HadamardProduct");
		{
			auto v = Create<ColVector<float>>(3, 1);
			auto m = Create<Matrix<float>>(3, 3);

			m = { {1,2,3}, {0, 3,4}, {2,3,4} };
			v = { 1,2,3 };
			

			std::cout << "m = " << m << std::endl;
			std::cout << "v = " << v << std::endl;

			std::cout << "ToMatrix(m) * ToMatrix(v) = " << ToMatrix(m) * ToMatrix(v) << std::endl;
			

			std::cout << "ScalarProduct(v, v) = "<< ScalarProduct(v, v) << std::endl;

			std::cout << "Transpose(m) = " << Transpose(m) << std::endl;

			std::cout << "HadamardProduct(m, m) = " << HadamardProduct(m, m) << std::endl;
		}

		counter("RowReduce, Rg, Det");
		{
			Matrix<double> m =
			{ {1., 0., 0.},
			  {4., 0., 0.},
			  {1., 0., 0.} 
			};
		

			double det;
			size_t rank;
			std::cout << RowReduce(m, &det, &rank) << std::endl;
			std::cout << "Det = " << det << std::endl;
			std::cout << "Rg = " << rank << std::endl;

			std::cout << "Det = " << Det(m) << std::endl;
			std::cout << "Rg = " << Rg(m) << std::endl;

		}

		counter("EuclidNorm, MaxNorm, FrobeniusNorm");
		{
			RowVector<int> v = { 1,2,3 };
			std::cout << "v = " << v << std::endl;

			std::cout << "EuclidNorm = " << EuclidNorm(v) << std::endl;

			std::cout << "MaxNorm = " << MaxNorm(v) << std::endl;

			std::cout << "v^T * v = " << Transpose(v) * v << std::endl;

			std::cout << "FrobeniusNorm(v^T * v) = " << FrobeniusNorm(Transpose(v) *  v) << std::endl;
		}
		
		counter("Angle");
		{
			RowVector<int> a = { 1,0,1 };
			ColVector<int> b = { 0,0,1 };

			std::cout << a << std::endl;
			std::cout << b << std::endl;
			std::cout << "deg = " << Angle(a, b, DEG) << std::endl;
			std::cout << "rad = " << Angle(a, b, RAD) << std::endl;
		}

		counter("SubMatrix, Assign, Inverse");
		{
			Matrix<float> m =
			{
				{1, 0},
				{0, 1},
				{4, 1}
			};

			ColVector<float> v = { 1, 2, 3 };
			auto a = Assign(m, v);
			std::cout << a;

			auto b = Assign(SubMatrix(a, EntryPosition{ 1, 0 }, EntryPosition{ 2, 1 }), ColVector<float>{ 1,2 });
			std::cout << b;

			Matrix<float> inv = Inverse(a); 
			std::cout <<"m^(-1) = " << inv << std::endl;


			std::cout <<"m * m^(-1) = " << a * inv << std::endl;
		}
	}
	end_section();


	begin_section("Дополнительный тест.");
	{
		Matrix<double> m =
		{ {1., 0., 0.},
		  {4., 5., 6.},
		};

		std::cout << Transpose(m);
		

		RowVector<int> rv = { 1,2,3 };
		std::cout << Transpose(rv);

	}
	end_section();


	//=====================================================
}
