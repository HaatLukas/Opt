#include "opt_alg.h"
#include <random>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>

// Funkcja do losowania wartoœci z zakresu [-2,2]
double randomForce(double minF, double maxF) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(minF, maxF);
	return dis(gen);
}

matrix fitnessFunction(matrix forces, matrix ud1, matrix ud2) {
	double m = 1.0;
	double b = 0.5;
	double D = 5.0;
	double T = 5.0;
	int n = get_len(forces);
	double dt = 0.01;
	int steps = static_cast<int>(T / dt);

	double x = 0.0, v = 0.0;
	double t = 0.0;
	int forceIndex = 0;
	double forceTime = T / n;

	double penalty = 0.0;

	for (int i = 0; i < steps; ++i) {
		if (t >= (forceIndex + 1) * forceTime && forceIndex < n - 1) {
			forceIndex++;
		}

		double F = forces(forceIndex, 0);

		if (F < -2 || F > 2) {
			penalty += 50.0 * (std::abs(F) - 2);  
		}

		double a = (F - b * v) / m;
		v += a * dt;
		x += v * dt;
		t += dt;
	}

	double error_x = std::abs(x - D);
	double error_v = std::abs(v);

	if (error_v > 0.3) {
		penalty += 3.0 * error_v;
	}

	matrix error(1, 1);
	error(0, 0) = error_x + 3 * error_v + penalty;

	std::cout << "Obliczanie funkcji celu dla F: " << forces << " -> Wynik: " << error << std::endl;

	return error;
}


// Algorytm genetyczny
solution GA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int pop_size, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	try {
		// Inicjalizacja populacji
		std::vector<solution> population;
		for (int i = 0; i < pop_size; ++i) {
			matrix forces(N, 1);
			for (int j = 0; j < N; ++j) {
				forces(j, 0) = randomForce(lb(j, 0), ub(j, 0));
			}
			solution sol(forces);
			sol.fit_fun(ff, ud1, ud2);
			population.push_back(sol);
		}

		int iter = 0;
		while (iter < Nmax) {
			// Selekcja metod¹ ko³a ruletki
			std::vector<solution> selected_population;
			double total_fitness = 0.0;
			for (const auto& sol : population) total_fitness += 1.0 / (sol.y(0, 0) + 1e-6);

			for (int i = 0; i < pop_size; ++i) {
				double rand_val = randomForce(0, total_fitness);
				double sum = 0.0;
				for (const auto& sol : population) {
					sum += 1.0 / (sol.y(0, 0) + 1e-6);
					if (sum >= rand_val) {
						selected_population.push_back(sol);
						break;
					}
				}
			}

			// Krzy¿owanie jednopunktowe
			for (int i = 0; i < pop_size; i += 2) {
				if (i + 1 < pop_size) {
					int crossover_point = rand() % N;
					for (int j = crossover_point; j < N; ++j) {
						std::swap(selected_population[i].x(j, 0), selected_population[i + 1].x(j, 0));
					}
				}
			}

			// Mutacja
			for (auto& sol : selected_population) {
				if (randomForce(0, 1) < 0.05) { // 5% szans na mutacjê
					int mut_index = rand() % N;
					double mutationAmount = randomForce(-0.3, 0.3);  // Mniejsza zmiana
					sol.x(mut_index, 0) = std::clamp(sol.x(mut_index, 0) + mutationAmount, -2.0, 2.0);
				}
			}

			for (auto& sol : selected_population) sol.fit_fun(ff, ud1, ud2);
			population = selected_population;

			iter++;
		}

		solution best_sol = population[0];
		for (const auto& sol : population) {
			if (sol.y(0, 0) < best_sol.y(0, 0)) best_sol = sol;
		}
		return best_sol;
	}
	catch (std::string ex_info) {
		throw ("solution GA(...):\n" + ex_info);
	}
}

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		//Tu wpisz kod funkcji

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
