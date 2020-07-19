#include<iostream>
#include<string>
#include<fstream>
#include<math.h>
using namespace std;

typedef struct Zuobiao
{
	double x;
	double y;
}Zuobiao;

double *duwenjian()
{
	int Num;//坐标总数
	string m;
	int i = 0, j;
	cout << "请输入数据文件名(以及路径)" << endl;
	cin >> m;
	fstream in;
	in.open(m, ios::in);
	do
	{
		if (i == 3)
		{
			in >> m;
			in >> Num;
		}
		getline(in, m);
		i++;
	} while (m != "NODE_COORD_SECTION");

	Zuobiao *matrix = new Zuobiao[Num + 1];
	double *Adjacency_matrix = new double[Num*(Num - 1) / 2 + 1];

	i = 1;
	while (i <= Num)
	{
		in >> m;
		if (m == "EOF")
			break;
		in >> matrix[i].x >> matrix[i].y;
		i++;
	}
	double temp;
	Adjacency_matrix[0] = Num;
	for (j = 1; j < Num + 1; j++)
	{
		for (i = j + 1; i < Num + 1; i++)
		{
			temp = sqrt((matrix[j].x - matrix[i].x)*(matrix[j].x - matrix[i].x) + (matrix[j].y - matrix[i].y)*(matrix[j].y - matrix[i].y));
			Adjacency_matrix[i*(i - 3) / 2 + j + 1] = temp;
		}
	}
	return Adjacency_matrix;
}
