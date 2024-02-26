#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
void DoOpen(double**& slau, double**& mA, int& size);
void DoPrint(double** slau, int size);
void Gauss(double**& slau, int& p, int size);
void ObrHod(double*& xx, double** slau, int size);
void FindNevyaz(double*& nevyaz, double** mA, double* xx, int size);
void ObrMatrx(double**& slau, int size);
void DoOut(double** slau, double* xx, double* nevyaz, int size, double opr);

int main()
{
	setlocale(LC_ALL, "RUSSIAN");
	int size;//размерность матрицы
	double** slau;//матрица имеющая вид (A)(b)(E)
	double** mA;//матрица А для невязок
	int perest = 0;//кол-во перестановок
	double* xx;//корни СЛАУ
	double* nevyaz; //невязки
	DoOpen(slau, mA, size); //функция чтения из файла
	DoPrint(slau, size); //Вывод матрицы
	xx = new double[size];
	Gauss(slau, perest, size); //Метод Гаусса
	ObrHod(xx, slau, size); //обратный ход
	FindNevyaz(nevyaz, mA, xx, size);//находим невязки
	double opr = pow(-1, perest); // находим определитель
	for (int i = 0; i < size; i++)
	{
		opr *= slau[i][i];
	}
	ObrMatrx(slau, size); //обратная матрица
	DoOut(slau, xx, nevyaz, size, opr); //Выведем матрицу имеющую вид (треугольный вид)(b)(Обратная матрица)
	system("pause");
	return 0;
}

void DoOpen(double**& slau, double**& mA, int& size)
{
	//Создаем файловый поток и связываем его с файлом
	ifstream file_input("input.txt");
	//Если файл открыт
	if (file_input.is_open())
	{
		file_input >> size;//считываем размерность квадратной матрицы А
		slau = new double* [2 * size + 1];//Размерность матрицы имеющей вид (A)(b)(E)
		mA = new double* [2 * size + 1];//Размерность матрицы А для невязок
		for (int i = 0; i < size; i++)
		{
			slau[i] = new double[2 * size + 1];
			mA[i] = new double[size + 1];
		}
		//Считаем матрицу из файла
		//Матрица А
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				file_input >> slau[i][j];
				mA[i][j] = slau[i][j];
			}
		//Вектор B
		for (int i = 0; i < size; i++)
		{
			file_input >> slau[i][size];
			mA[i][size] = slau[i][size];
		}
		file_input.close();//закрываем файл
	}
	else
	{
		//Если открытие файла прошло не успешно
		cout << "Файл не найден.";
		exit(-1);
	}
	//строим единичную матрицу
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i == j)
				slau[i][size + 1 + j] = 1;
			else
				slau[i][size + 1 + j] = 0;
		}
	}
}
void DoPrint(double** slau, int size)
{
	cout << "Матрица имеющая вид (A)(b)(E)" << endl;
	//Выведем матрицу имеющую вид (A)(b)(E)
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 2 * size + 1; j++)
			cout << slau[i][j] << "\t ";
		cout << endl;
	}
}
void Gauss(double**& slau, int& p, int size)
{
	double max;//макс элемент в столбце
	int index;//индекс макс элемента
	//Прямой ход, приведение к верхнетреугольному виду
	double tmp;//число, на которое делим
	const double eps = 0.0001;// точность

	for (int i = 0; i < size; i++)
	{
		max = abs(slau[i][i]);
		index = i;
		//поиск макс элемента
		for (int k = i + 1; k < size; k++)
		{
			if (abs(slau[k][i]) > max)
			{
				max = abs(slau[k][i]);
				index = k;
			}
		}
		//Определение вырожденности матрицы
		if (max < eps)
		{
			cout << "Решение получить невозможно, матрица вырожденная!";
			exit(0);
		}
		//перестановка строк
		if (index != i)
		{
			for (int j = 0; j < 2 * size + 1; j++)
			{
				double buff = slau[i][j];
				slau[i][j] = slau[index][j];
				slau[index][j] = buff;
			}
			p++;
		}
		//прямой ход
		for (int k = i + 1; k < size; k++) {
			tmp = -slau[k][i] / slau[i][i];
			for (int j = 0; j < size + 1; j++) {
				slau[k][j] += slau[i][j] * tmp;
				slau[k][j + size + 1] += slau[i][j + size + 1] * tmp;
			}
		}
	}
}
void ObrHod(double*& xx, double** slau, int size)
{
	xx[size - 1] = slau[size - 1][size] / slau[size - 1][size - 1];
	for (int i = size - 2; i >= 0; i--)
	{
		xx[i] = slau[i][size];
		for (int j = i + 1; j < size; j++)
			xx[i] -= slau[i][j] * xx[j];
		xx[i] /= slau[i][i];
	}
}
void FindNevyaz(double*& nevyaz, double** mA, double* xx, int size)
{
	nevyaz = new double[size];
	for (int i = 0; i < size; i++) {
		nevyaz[i] = mA[i][size];
		for (int j = 0; j < size; j++) {
			nevyaz[i] -= mA[i][j] * xx[j];
		}
	}
	cout << endl << "Невязки" << endl;
	for (int i = 0; i < size; i++) {
		cout << nevyaz[i] << " ";
	}
	cout << endl;
	
}
void ObrMatrx(double**& slau, int size)
{
	for (int k = 0; k < size; k++) {
		slau[size - 1][size + 1 + k] = slau[size - 1][size + 1 + k] / slau[size - 1][size - 1];
		for (int i = size - 2; i >= 0; i--)
		{
			slau[i][size + 1 + k] = slau[i][size + 1 + k];
			for (int j = i + 1; j < size; j++)
				slau[i][size + 1 + k] -= slau[i][j] * slau[j][size + 1 + k];
			slau[i][size + 1 + k] /= slau[i][i];
		}
	}
}
void DoOut(double** slau, double* xx, double* nevyaz, int size, double opr)
{
	cout << endl << "Матрица имеющая вид (треугольная вид)(b)(Обратная матрица)" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 2 * size + 1; j++)
			cout << slau[i][j] << "\t ";
		cout << endl;
	}
	cout << endl << "Решение системы" << endl;
	for (int i = 0; i < size; i++)
		cout << xx[i] << " ";
	cout << endl;
	cout << endl << "Определитель" << endl;
	cout << opr;
	cout << endl;
	ofstream file_output;
	file_output.open("output.txt");
	if (file_output.is_open())
	{
		for (int i = 0; i < size; i++)
			file_output << xx[i] << "\t ";
		file_output << endl;
		for (int i = 0; i < size; i++) {
			file_output << nevyaz[i] << "\t ";
		}
		file_output << endl;
		for (int k = 0; k < size; k++) {
			for (int g = size + 1; g < size * 2 + 1; g++) {
				file_output << slau[k][g] << "\t ";
			}
			file_output << endl;
		}
		file_output << endl;
		file_output << opr << endl;;
		file_output.close();//закрываем файл
	}
	else
	{
		//Если открытие файла прошло не успешно
		cout << "Файл не найден.";
		exit(-1);
	}
}

