#include<iostream>
#include<cstdlib>
#include<vector>
#include<string>
#include<time.h>
#include<fstream>
#include<math.h>
#include<sys/timeb.h>
using namespace std;
double* Adjacency_matrix;//压缩后的邻接矩阵
typedef struct Zuobiao
{
	double x;
	double y;
}Zuobiao;
double get_distance(int i, int j)
{
	if (i > j)
		return Adjacency_matrix[i*(i - 3) / 2 + j + 1];
	else
		return Adjacency_matrix[j*(j - 3) / 2 + i + 1];
}
__declspec (naked) unsigned __int64 GetCpuCycle(void)
{
	_asm
	{
		rdtsc
		ret
	}
}
int get_rand(int x)//生成1-x(含)以内的随机数
{
	unsigned __int64 iCpuCycle = GetCpuCycle();
	unsigned srnd = (unsigned)iCpuCycle;
	srand(srnd);
	return (rand() % x + 1);
}

void duwenjian(string m)
{
	int Num;//坐标总数
	int i = 0, j;
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
	} 
	while (m != "NODE_COORD_SECTION");
	Zuobiao *matrix = new Zuobiao[Num + 1];
	Adjacency_matrix = new double[Num*(Num - 1) / 2 + 1];
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
	return;
}

class GA
{
private:
	vector<vector<int>>fathers_gene_pool;
	vector<vector<int>>sons_gene_pool;
	vector<int>best_path;//当前代最优解
	vector<int>fitness;//适应度
	vector<int>arrival;//记录已经访问过的城市，访问过对应地点为1，否则为0
	int generation;//当前遗传代数
	int target;//目标遗传代数
	double P_mutation;//变异概率
	double P_cross;//交叉概率
	int population;//种群数量
public:
	GA(int t, int num, int population, double Pm, double Pc);//类构造函数
	void create(int population);//初始解的构造函数
	void mutation(int gene);//基因池中的第几个发生变异
	void cross(int gene1,int gene2);//基因池中的两个发生交叉
	double cal_fitness(int gene);//计算指定基因的适应度
	bool select();//根据遗传代数概率性接受较差的解
	void print(int type, int num);
	int* draw_path();
};
GA::GA(int t, int num, int p, double Pm, double Pc)//目标遗传代数 城市数量 初始种群大小
{
	target = t;
	vector<int> temp(num + 1, 0);
	temp[1] = 1;
	arrival = temp;
	P_cross = Pc;
	P_mutation = Pm;
	population = p;
	generation = 0;
	create(p);
}
void GA::create(int population)
{
	vector<int>temp;
	for (int i = 0; i < population; i++)
	{
		fathers_gene_pool.push_back(temp);
		fathers_gene_pool[i].push_back(1);
	}
	for (int i = 1; i <= Adjacency_matrix[0]; i++)
	{
		double nearest_d = INFINITY;
		int nearest_p = 0;
		for (int j = 2; j <= Adjacency_matrix[0]; j++)
		{
			if (arrival[j])
				continue;
			else
			{
				double d = get_distance(i, j);
				if (nearest_d > d)
				{
					nearest_d = d;
					nearest_p = j;
				}
			}
		}
		if (nearest_p == 0)
			fathers_gene_pool[0].push_back(1);
		else
		{
			fathers_gene_pool[0].push_back(nearest_p);
			arrival[nearest_p] = 1;
		}
	}
	for (int i = 1; i < population; i++)
	{
		fathers_gene_pool[i].assign(fathers_gene_pool[0].begin(), fathers_gene_pool[0].end());
	}

	for (int j = 1; j < population; j++)
	{
		for (int i = 0; i < 10; i++)
		{
			int pos1 = get_rand(fathers_gene_pool[0].size() - 2);
			int pos2 = get_rand(fathers_gene_pool[0].size() - 2);
			swap(fathers_gene_pool[j][pos1], fathers_gene_pool[j][pos2]);
		}
	}
}
void GA::mutation(int gene)
{
	int pos1 = get_rand(Adjacency_matrix[0]-1);
	int pos2 = get_rand(Adjacency_matrix[0]-1);
	while(pos1 == pos2)
		pos2 = get_rand(Adjacency_matrix[0]-1);
	swap(fathers_gene_pool[gene][pos1], fathers_gene_pool[gene][pos2]);
	sons_gene_pool.push_back(fathers_gene_pool[gene]);
}
void GA::cross(int gene1, int gene2)
{
	int pos1 = get_rand(Adjacency_matrix[0] - 1);
	int pos2 = get_rand(Adjacency_matrix[0] - 1);
	//print(1, gene1);
	//print(1, gene2);
	//cout << endl;
	int start = pos1 < pos2 ? pos1 : pos2;
	int end = start + abs(pos1 - pos2);
	for (int i = start; i <= end; i++)
		swap(fathers_gene_pool[gene1][i], fathers_gene_pool[gene2][i]);
	/*print(1, gene1);
	print(1, gene2);
	cout << endl;*/
	int times = 0;
	for (int i = start; i <= end; i++)
	{
		vector<int>conflict1;
		vector<int>conflict2;
		for (int i = start; i <= end; i++)
		{
			times = count(fathers_gene_pool[gene1].begin(), fathers_gene_pool[gene1].end(), fathers_gene_pool[gene1][i]);
			if (times == 2)
				for (int j = 0; j < fathers_gene_pool[gene1].size(); j++)
					if ((fathers_gene_pool[gene1][j] == fathers_gene_pool[gene1][i]) && (j != i))
						conflict1.push_back(j);
		}
		for (int i = start; i <= end; i++)
		{
			times = count(fathers_gene_pool[gene2].begin(), fathers_gene_pool[gene2].end(), fathers_gene_pool[gene2][i]);
			if (times == 2)
				for (int j = 0; j < fathers_gene_pool[gene2].size(); j++)
					if ((fathers_gene_pool[gene2][j] == fathers_gene_pool[gene2][i]) && (j != i))
						conflict2.push_back(j);
		}
		for (int i = 0; i < conflict1.size(); i++)
		{
			swap(fathers_gene_pool[gene1][conflict1[i]], fathers_gene_pool[gene2][conflict2[i]]);
		}
	}
	/*print(1, gene1);
	print(1, gene2);
	cout << endl;*/
	sons_gene_pool.push_back(fathers_gene_pool[gene1]);
	sons_gene_pool.push_back(fathers_gene_pool[gene2]);
}
double GA::cal_fitness(int gene)
{
	double fitness = 0;
	for (int i = 0; i < sons_gene_pool[gene].size() - 2; i++)
	{
		fitness += get_distance(sons_gene_pool[gene][i], sons_gene_pool[gene][i + 1]);
	}
	return fitness;
}
void GA::print(int type, int num)//1fu
{
	if(type == 1)
		for (int i = 0; i < fathers_gene_pool[num].size(); i++)
			cout << fathers_gene_pool[num][i]<<" ";
	else
		for (int i = 0; i < sons_gene_pool[num].size(); i++)
			cout << sons_gene_pool[num][i]<<" ";
	cout << endl;
}
int* GA::draw_path()
{
	while (generation < target)
	{
		sons_gene_pool.clear();//清空子代准备新一轮
		sons_gene_pool.push_back(fathers_gene_pool[fathers_gene_pool.size() - 1]);
		for (int i = 0; i < population; i++)//遗传变异
		{
			if (get_rand(1000) < P_cross * 1000) 
			{
				int t = get_rand(population) - 1;
				cross(i,t);
			}
			if (get_rand(1000) < P_mutation * 1000)
				mutation(i);
			if (get_rand(1000) < P_mutation * 1000)
				mutation(i);
		}
		vector<double>fit;
		for (int i = 0; i < sons_gene_pool.size(); i++)//计算适应度
		{
			fit.push_back(cal_fitness(i));
		}
		fathers_gene_pool.clear();//清空父本准备新一轮
		double min = INFINITY;
		int min_pos = 0;
		while (fathers_gene_pool.size() < population)//筛选质量好的子代，二元锦标赛法
		{
			int pos1 = get_rand(sons_gene_pool.size()) - 1;
			int pos2 = get_rand(sons_gene_pool.size()) - 1;
			if (fit[pos1] > fit[pos2])
			{
				vector<int>temp = sons_gene_pool[pos2];
				fathers_gene_pool.push_back(temp);
			}
			else
				fathers_gene_pool.push_back(sons_gene_pool[pos1]);
			if (fit[pos2] < min)
			{
				min = fit[pos2];
				min_pos = pos2;
			}
			else if (fit[pos1] < min)
			{
				min = fit[pos1];
				min_pos = pos1;
			}
		}
		generation++;
		cout << generation << endl;
		cout <<"当前代最优解"<< min << endl;
		fathers_gene_pool.push_back(sons_gene_pool[min_pos]);
		print(0, min_pos);
	}
	return nullptr;
}
int main(void)
{
	string filename;
	cout << "请输入文件名" << endl;
	cin >> filename;
	filename += ".tsp";
	duwenjian(filename);
	GA ga(1000, Adjacency_matrix[0],100,0.1,0.8);
	ga.draw_path();
	cout << endl;
	return 0;
}