#include <iostream>
#include <windows.h>  // For SetConsoleOutputCP
#include "TrajectoryOptimizer.hpp"

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    // Создаем оптимизатор траектории
    TrajectoryOptimizer optimizer(1.0, 0.1, 20);
    
    std::cout << "=== Оптимизация траекторий ===" << std::endl;
    
    // Оптимизация по расходу топлива
    std::cout << "\nОптимизация по расходу топлива:" << std::endl;
    DMDict result_fuel = optimizer.trajectory("fuel", {{"xf", 5.0}, {"yf", 3.0}});
    optimizer.visualizeTrajectory(result_fuel);
    
    // Оптимизация по энергии
    std::cout << "\nОптимизация по энергии:" << std::endl;
    DMDict result_energy = optimizer.trajectory("energy", {{"xf", 5.0}, {"yf", 3.0}});
    optimizer.visualizeTrajectory(result_energy);
    
    return 0;
}