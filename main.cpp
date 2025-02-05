#include "opt_alg.h"
#include <iostream>

int main() {
    int pop_size = 50;
    double epsilon = 1e-3;
    int Nmax = 10;

    for (int n = 2; n <= 5; ++n) 
    {
        std::cout << "Testowanie dla n = " << n << ":\n";

        matrix lb(n, 1, -2.0);  // Teraz lb ma wymiar n x 1
        matrix ub(n, 1, 2.0);   // Teraz ub ma wymiar n x 1

        for (int i = 0; i < 1; ++i) {
            try {
                solution result = GA(fitnessFunction, n, lb, ub, pop_size, epsilon, Nmax, NAN, NAN);
                std::cout << "Iteracja " << i + 1 << ": y* = " << result.y(0, 0) << "\n";
            }
            catch (const std::string& ex) {
                std::cerr << "B³¹d: " << ex << std::endl;
            }
            catch (const std::exception& ex) {
                std::cerr << "Wyj¹tek: " << ex.what() << std::endl;
            }
            catch (...) {
                std::cerr << "Nieznany b³¹d!" << std::endl;
            }
        }
        std::cout << "--------------------------------\n";
    }
    return 0;
}
