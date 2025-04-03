#include <iostream>
#include <windows.h>  // For SetConsoleOutputCP
#include <memory>
#include "TrajectoryOptimizer.hpp"


void demoAPI() {
    std::cout << "\n=== ДЕМОНСТРАЦИЯ API ===" << std::endl;
    
    // Создаем параметры системы - объект массой 2 кг с ограничением силы 3 Н
    SystemParams system_params(2.0, 3.0);
    
    // Создаем параметры оптимизации - 30 шагов по 0.2 сек с использованием IPOPT
    OptimizationParams opt_params(30, 0.2, "ipopt");
    
    // Создаем оптимизатор с нашими параметрами
    std::unique_ptr<TrajectoryOptimizer> optimizer = 
        std::make_unique<TrajectoryOptimizer>(system_params, opt_params);
    
    // Оптимизация по расходу топлива
    std::cout << "\nОптимизация по расходу топлива:" << std::endl;
    TrajectoryResult result_fuel = optimizer->optimize("fuel", {{"xf", 8.0}, {"yf", 5.0}});
    optimizer->visualizeTrajectory(result_fuel);
    
    // Изменяем параметры системы - теперь объект легче
    SystemParams light_system(0.5, 3.0);
    optimizer->setSystemParams(light_system);
    
    // Оптимизация по энергии для легкого объекта
    std::cout << "\nОптимизация по энергии для легкого объекта (масса = 0.5 кг):" << std::endl;
    TrajectoryResult result_energy_light = optimizer->optimize("energy", {{"xf", 8.0}, {"yf", 5.0}});
    optimizer->visualizeTrajectory(result_energy_light);
    
    // Создаем другой оптимизатор для демонстрации минимизации времени
    std::unique_ptr<TrajectoryOptimizer> time_optimizer = std::make_unique<TrajectoryOptimizer>();
    
    // Настраиваем параметры
    time_optimizer->setSystemParams(SystemParams(1.0, 4.0));
    time_optimizer->setOptimizationParams(OptimizationParams(40, 0.05, "ipopt"));
    
    // Оптимизация по времени
    std::cout << "\nОптимизация по времени:" << std::endl;
    TrajectoryResult result_time = time_optimizer->optimize("time", {{"xf", 10.0}, {"yf", 10.0}});
    time_optimizer->visualizeTrajectory(result_time);
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    demoAPI();

    return 0;
}