#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <stdio.h>
#include <process.h>
#include <stdlib.h>
#include <Windows.h>
#include <ctime>
#include <deque>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <iostream>

#define sqr(x) ((x)*(x))

using namespace std;

struct people
{
	int xs, ys, xf, yf, ind, time;
	people(int xs, int ys, int xf, int yf, int ind, int time) : xs(xs), ys(ys), xf(xf), yf(yf), ind(ind), time(time){};
};

struct people_data
{
	long long int T_in, T_out, T_start, T_best;
	people P;
	people_data() : T_in(0), T_out(0), T_start(0), T_best(0), P(people(0, 0, 0, 0, 0, 0)){}
};

struct step
{
	int x, y, type;
	step(int new_x, int new_y, int type) : x(new_x), y(new_y), type(type){}
};

struct TAXI
{
	int x, y, free_places;
	TAXI(int new_x, int new_y, int free_places) : x(new_x), y(new_y), free_places(free_places){}
	TAXI() : x(0), y(0), free_places(4){}
};



int delta_to_make_everyone_unhappy_happy = 500;
int size_coef_to_make_everyone_unhappy_happy = 10;
double delta_border = 0.5;
int eps = 150;
double eps_coef = 2;
double additional_passenger_coef = 0.5;
int find_companion_border = 500;

int mid_sum = 0;
int mid_dist = 0;

bool used[1000];
vector<TAXI> taxi;
vector<pair<int, int>> good_pos;
deque <step> go[1000], additional_go[1000];
people_data p_d[1000];
int penalty[1000];
int free_taxi_T[1000];

vector<people> P;
int w, h, n;
int global_time = 0;

void clear_q(deque<step> &Q)
{
	while (!Q.empty())
		Q.pop_back();
}

void kik_passenger(int index)
{
	for (int i = 0; i < P.size(); i++)
		if (P[i].ind == index)
			P.erase(P.begin() + i);
}



double calc_profit_main(int ind)
{
	return 0;
	TAXI T = taxi[ind];
	deque<step> GO = go[ind];
	int dist = 0;
	double res = 0;
	int new_T = 0;
	for (int i = 0; i < GO.size(); i++)
	{
		if (i == 0)
			dist += abs(T.x - GO[i].x) + abs(T.y + GO[i].y);
		else
			dist += abs(GO[i - 1].x - GO[i].x) + abs(GO[i - 1].y + GO[i].y);
		new_T = global_time + dist;
		if (GO[i].type > 0)
			p_d[GO[i].type - 1].T_in = new_T;
		else
		{
			res += (1 - (sqr(p_d[GO[i].type - 1].T_in - p_d[GO[i].type - 1].T_start + 0.0) +
				sqr(p_d[i].T_best - new_T + p_d[i].T_in)) / 10000000.0) * 
				(abs(p_d[GO[i].type - 1].P.xs - p_d[GO[i].type - 1].P.xf) + abs(p_d[GO[i].type - 1].P.ys - p_d[GO[i].type - 1].P.yf));
		}
	}
	return res;
}


double calc_profit_all(int ind)
{
	return 0;
	TAXI T = taxi[ind];
	deque<step> GO = go[ind];
	for (int i = 0; i < additional_go[ind].size(); i++)
	{
		GO.push_back(additional_go[ind].front());
		additional_go[ind].push_back(additional_go[ind].front());
		additional_go[ind].pop_front();
	}
	int dist = 0;
	double res = 0;
	int new_T = 0;
	for (int i = 0; i < GO.size(); i++)
	{
		if (i == 0)
			dist += abs(T.x - GO[i].x) + abs(T.y + GO[i].y);
		else
			dist += abs(GO[i - 1].x - GO[i].x) + abs(GO[i - 1].y + GO[i].y);
		new_T = global_time + dist;
		if (GO[i].type > 0)
			p_d[GO[i].type - 1].T_in = new_T;
		else
		{
			res += (1 - (sqr(p_d[GO[i].type - 1].T_in - p_d[GO[i].type - 1].T_start + 0.0) +
				sqr(p_d[i].T_best - new_T + p_d[i].T_in)) / 10000000.0) *
				(abs(p_d[GO[i].type - 1].P.xs - p_d[GO[i].type - 1].P.xf) + abs(p_d[GO[i].type - 1].P.ys - p_d[GO[i].type - 1].P.yf));
		}
	}
	return res;
}


void print_all()
{
	int cnt = 0;
	cout << n << ' ';
	for (int i = 0; i < n; i++)
	{
		cout << i + 1 << " " << go[i].size() + additional_go[i].size() << " ";
		for (int j = 0; j < go[i].size(); j++)
		{
			cout << go[i].front().x << ' ' << go[i].front().y << ' ' << go[i].front().type << ' ';
			go[i].push_back(go[i].front());
			go[i].pop_front();
		}
		for (int j = 0; j < additional_go[i].size(); j++)
		{
			cout << additional_go[i].front().x << ' ' << additional_go[i].front().y << ' ' << additional_go[i].front().type << ' ';
			additional_go[i].push_back(additional_go[i].front());
			additional_go[i].pop_front();
		}
	}
	cout << endl;
}

void simulate(int delta_t)
{
	for (int i = 0; i < n; i++)
	{
		int T = delta_t;
		while (!go[i].empty())
		{
			if (T >= abs(taxi[i].x - go[i].front().x) + abs(go[i].front().y - taxi[i].y))
			{
				while (!go[i].empty() && T >= abs(taxi[i].x - go[i].front().x) + abs(go[i].front().y - taxi[i].y))
				{
					T -= abs(taxi[i].x - go[i].front().x) + abs(go[i].front().y - taxi[i].y);
					taxi[i].x = go[i].front().x;
					taxi[i].y = go[i].front().y;
					if (go[i].front().type < 0)
					{
						taxi[i].free_places++;
						p_d[abs(go[i].front().type) - 1].T_out = global_time - T + delta_t;
					}
					else
						if (go[i].front().type > 0)
						{
							p_d[abs(go[i].front().type) - 1].T_in = global_time - T + delta_t;
						}
					go[i].pop_front();
				}
			}
			else
				if (T >= abs(go[i].front().x - taxi[i].x))
				{
					T -= abs(taxi[i].x - go[i].front().x);
					taxi[i].x = go[i].front().x;

					int dir = (go[i].front().y - taxi[i].y) / abs(go[i].front().y - taxi[i].y);
					taxi[i].y += T * dir;
					T = 0;
				}
				else
				{
					int dir = (go[i].front().x - taxi[i].x) / abs(go[i].front().x - taxi[i].x);
					taxi[i].x += T * dir;
					T = 0;
				}
			if (!T)
				break;
		}
		if (go[i].empty())
		{
			while (!additional_go[i].empty())
			{
				if (T >= abs(taxi[i].x - additional_go[i].front().x) + abs(additional_go[i].front().y - taxi[i].y))
				{
					while (!additional_go[i].empty() && T >= abs(taxi[i].x - additional_go[i].front().x) + abs(additional_go[i].front().y - taxi[i].y))
					{
						T -= abs(taxi[i].x - additional_go[i].front().x) + abs(additional_go[i].front().y - taxi[i].y);
						taxi[i].x = additional_go[i].front().x;
						taxi[i].y = additional_go[i].front().y;
						if (additional_go[i].front().type < 0)
							p_d[abs(additional_go[i].front().type) - 1].T_out = global_time - T + delta_t;

						if (additional_go[i].front().type > 0)
						{
							p_d[abs(additional_go[i].front().type) - 1].T_in = global_time - T + delta_t;
							kik_passenger(additional_go[i].front().type);
						}
						additional_go[i].pop_front();
					}
				}
				else
					if (T >= abs(additional_go[i].front().x - taxi[i].x))
					{
						T -= abs(taxi[i].x - additional_go[i].front().x);
						taxi[i].x = additional_go[i].front().x;
						int dir = (additional_go[i].front().y - taxi[i].y) / abs(additional_go[i].front().y - taxi[i].y);
						taxi[i].y += T * dir;
						T = 0;
					}
					else
					{
						int dir = (additional_go[i].front().x - taxi[i].x) / abs(additional_go[i].front().x - taxi[i].x);
						taxi[i].x += T * dir;
						T = 0;
					}
				if (!T)
					break;
			}

			if (!additional_go[i].empty())
				if (additional_go[i].front().type < 0)
				{
					go[i].push_back(additional_go[i].front());
					additional_go[i].pop_front();
					taxi[i].free_places--;
				}
		}
	}
}

void make_everyone_unhappy_happy()
{
	double result_T[1000];

	for (int j = 0; j < n; j++)
	{
		additional_go[j].clear();
		result_T[j] = calc_profit_main(j);
	}

	for (int i = P.size() - 1; i >= 0; i--)
	{

		int best = LONG_MAX;
		int ind = -1;
		for (int j = 0; j < n; j++)
		{
			double delta = (1 - sqr(global_time - P[i].time + 0.0) / 10000000.0);
			additional_go[j].push_back(step(P[i].xs, P[i].ys, P[i].ind));
			additional_go[j].push_back(step(P[i].xf, P[i].yf, -P[i].ind));
			if (delta > delta_border)
			{

				if (go[j].size() == 0)
				{
					if (abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) < best)
					{
						ind = j;
						best = abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys);
					}
				}
				else
					if (additional_go[j].size() == 0)
					{
						if (abs(go[j].back().x - P[i].xs) + abs(go[j].back().y - P[i].ys) + go[j].size() * (w + h) / size_coef_to_make_everyone_unhappy_happy < best)
						{
							ind = j;
							best = abs(go[j].back().x - P[i].xs) + abs(go[j].back().y - P[i].ys) + go[j].size() * (w + h) / size_coef_to_make_everyone_unhappy_happy;
						}
					}
					else
						if (abs(additional_go[j].back().x - P[i].xs) + abs(additional_go[j].back().y - P[i].ys) + (go[j].size() + additional_go[j].size()) * (w + h) / size_coef_to_make_everyone_unhappy_happy < best)
						{
							ind = j;
							best = abs(additional_go[j].back().x - P[i].xs) + abs(additional_go[j].back().y - P[i].ys) + (go[j].size() + additional_go[j].size()) * (w + h) / size_coef_to_make_everyone_unhappy_happy;
						}
			}
			additional_go[j].pop_back();
			additional_go[j].pop_back();
		}
		if (ind != -1)
		{
			additional_go[ind].push_back(step(P[i].xs, P[i].ys, P[i].ind));
			additional_go[ind].push_back(step(P[i].xf, P[i].yf, -P[i].ind));
			result_T[ind] = calc_profit_all(ind);
		}
	}
}

void find_companion(int j)
{
	for (int i = P.size() - 1; i >= 0; i--)
	{
		int best = LONG_MAX;
		double delta = (1 - sqr(global_time - P[i].time + 0.0) / 10000000.0);
		if (delta > delta_border)
		{
			if (taxi[j].free_places > 0 && go[j].size() > 0)
			{
				if (((taxi[j].x - go[j].front().x <= 0 && taxi[j].x - P[i].xs <= 0) || (taxi[j].x - go[j].front().x >= 0 && taxi[j].x - P[i].xs >= 0)) && abs(taxi[j].x - P[i].xs) < abs(taxi[j].x - go[j].front().x)
					&& ((taxi[j].y - go[j].front().y <= 0 && taxi[j].y - P[i].ys <= 0) || (taxi[j].y - go[j].front().y >= 0 && taxi[j].y - P[i].ys >= 0)) && abs(taxi[j].y - P[i].ys) < abs(taxi[j].y - go[j].front().y))
				{
					if (abs(go[j].back().x - P[i].xf) + abs(go[j].back().y - P[i].yf) < find_companion_border)
					{
						taxi[j].free_places--;
						go[j].push_back(step(P[i].xf, P[i].yf, -P[i].ind));
						go[j].push_front(step(P[i].xs, P[i].ys, P[i].ind));
						P.erase(P.begin() + i);
						break;
					}
				}
			}
		}
	}
}



void find_best()
{
	bool additional_passenger = false;
	for (int i = P.size() - 1; i >= 0; i--)
	{
		int best = LONG_MAX;
		int ind = -1;
		if (P.size() > i)
		{

				for (int j = 0; j < n; j++)
				{
					double delta = (1 - sqr(global_time - P[i].time + 0.0 + abs(P[i].xs - taxi[j].x) + abs(P[i].ys - taxi[j].y)) / 10000000.0);
					if (delta > 0.1)
						//или если долго свободен
						if (abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) < 3 * max(abs(P[i].xf - P[i].xs) + abs(P[i].yf - P[i].ys), 300) && abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) < (w + h) / 3)
							if (taxi[j].free_places > 0)
							{
								if (taxi[j].free_places == 4 && (abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) < best))
								{
									ind = j;
									best = abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys);
									additional_passenger = false;
								}
								else
								{
									if (taxi[j].free_places < 4 && abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) < eps && abs(go[j].back().x - P[i].xf) + abs(go[j].back().y - P[i].yf) < eps_coef * eps && best > (abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys)) * additional_passenger_coef)
									{
										ind = j;
										best = (abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys)) * additional_passenger_coef;
										additional_passenger = true;
									}
								}
							}
				}
		}
		if (ind != -1)
		{
			taxi[ind].free_places--;
			go[ind].push_back(step(P[i].xf, P[i].yf, -P[i].ind));
			go[ind].push_front(step(P[i].xs, P[i].ys, P[i].ind));
			P.erase(P.begin() + i);
			find_companion(ind);
		}
	}
}


void make_everyone_unhappy_happy_last()
{
	for (int j = 0; j < n; j++)
		additional_go[j].clear();
	for (int i = P.size() - 1; i >= 0; i--)
	{
		int best = LONG_MAX;
		int ind = -1;
		for (int j = 0; j < n; j++)
		{
			if (go[j].size() > 0)
			{
				if (abs(go[j].back().x - P[i].xs) + abs(go[j].back().y - P[i].ys) + go[j].size() * (w + h) / size_coef_to_make_everyone_unhappy_happy < best)
				{
					ind = j;
					best = abs(go[j].back().x - P[i].xs) + abs(go[j].back().y - P[i].ys) + go[j].size() * (w + h) / size_coef_to_make_everyone_unhappy_happy;
				}
			}
			else
			{
				if (abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) + go[j].size() * (w + h) / size_coef_to_make_everyone_unhappy_happy < best)
				{
					ind = j;
					best = abs(taxi[j].x - P[i].xs) + abs(taxi[j].y - P[i].ys) + go[j].size() * (w + h) / size_coef_to_make_everyone_unhappy_happy;
				}
			}
		}
		if (ind != -1)
		{
			go[ind].push_back(step(P[i].xs, P[i].ys, P[i].ind));
			go[ind].push_back(step(P[i].xf, P[i].yf, -P[i].ind));
			P.pop_back();
		}
	}
}

void print_taxi()
{
	cout << endl << global_time << endl;
	for (int i = 0; i < n; i++)
		cout << taxi[i].x << ' ' << taxi[i].y << endl;
	cout << "#######";
}


int main(int argc, char *argv[])
{


#ifndef ONLINE_JUDGE
	freopen(argv[4], "r", stdin);
	freopen(argv[5], "w", stdout);
#endif
	cin >> w >> h;
	cin >> n;
	taxi.assign(n, TAXI());
	for (int i = 0; i < n; i++)
		cin >> taxi[i].x >> taxi[i].y;

	int ind = 0;
	bool fl_take_all = false;
	int cnt = 0;
	int TT = 0;
	int delta_t = 0;

	while (global_time >= 0)
	{
		cnt++;
		print_all();
		//print_taxi();
		ind++;
		int t, xs, ys, xf, yf;
		cin >> t >> xs >> ys >> xf >> yf;
		mid_sum = abs(xs - xf) + abs(ys - yf);
		mid_dist = mid_sum / cnt;
		if (t != -1)
		{
			p_d[ind - 1].T_start = t;
			p_d[ind - 1].T_best = abs(xs - xf) + abs(ys - yf);
			simulate(t - global_time);
			P.push_back(people(xs, ys, xf, yf, ind, t));
			p_d[ind].P = P.back();
			TT = t;
		}
		delta_t = max(t - global_time, delta_t);
		global_time = t;
		find_best();
		for (int i = 0; i < n; i++)
			find_companion(i);
		if (t != -1 && delta_t > delta_to_make_everyone_unhappy_happy) // Перебрать 500
			make_everyone_unhappy_happy();
	}
	make_everyone_unhappy_happy_last();
	print_all();
	global_time = TT;
	simulate(10000000);

	//for (int i = 0; i < cnt - 1; i++)
	//	if (p_d[i].T_in != 0)
	//		sum += (100 + p_d[i].T_best) * (10000000.0 - min(sqr(p_d[i].T_start - p_d[i].T_in) + sqr(p_d[i].T_best - p_d[i].T_out + p_d[i].T_in), 10000000.0)) / 10000000.0;

	//cout << sum / (cnt - 1);
}