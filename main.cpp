#define NOMINMAX

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
    TrajectoryResult result_fuel = optimizer->optimize2D("fuel", {{"xf", 8.0}, {"yf", 5.0}});
    optimizer->visualizeTrajectory(result_fuel);
    
    // Изменяем параметры системы - теперь объект легче
    SystemParams light_system(0.5, 3.0);
    optimizer->setSystemParams(light_system);
    
    // Оптимизация по энергии для легкого объекта
    std::cout << "\nОптимизация по энергии для легкого объекта (масса = 0.5 кг):" << std::endl;
    TrajectoryResult result_energy_light = optimizer->optimize2D("energy", {{"xf", 8.0}, {"yf", 5.0}});
    optimizer->visualizeTrajectory(result_energy_light);
    
    // Создаем другой оптимизатор для демонстрации минимизации времени
    std::unique_ptr<TrajectoryOptimizer> time_optimizer = std::make_unique<TrajectoryOptimizer>();
    
    // Настраиваем параметры
    time_optimizer->setSystemParams(SystemParams(1.0, 4.0));
    time_optimizer->setOptimizationParams(OptimizationParams(40, 0.05, "ipopt"));
    
    // Оптимизация по времени
    std::cout << "\nОптимизация по времени:" << std::endl;
    TrajectoryResult result_time = time_optimizer->optimize2D("time", {{"xf", 10.0}, {"yf", 10.0}});
    time_optimizer->visualizeTrajectory(result_time);
}

void demo4DTrajectory() {
    std::cout << "\n=== ДЕМОНСТРАЦИЯ 4D-ОПТИМИЗАЦИИ ТРАЕКТОРИИ ===" << std::endl;
    
    // Создаем параметры системы - самолет массой 80 тонн
    SystemParams aircraft_params(80000.0, 200000.0, 100000.0, 250.0, 3000.0, 12000.0);
    
    // Параметры оптимизации - 50 точек с использованием IPOPT
    OptimizationParams opt_params(50, 0.0, "ipopt"); // dt=0 означает, что время тоже оптимизируется
    
    // Создаем оптимизатор с нашими параметрами
    std::unique_ptr<TrajectoryOptimizer> optimizer = 
        std::make_unique<TrajectoryOptimizer>(aircraft_params, opt_params);
    
    // Добавляем ограничения воздушного пространства (грозовые очаги, запретные зоны)
    
    // Грозовой фронт - сферическая область
    optimizer->addConstraint(AirspaceConstraint::createSphere(
        50000.0, 40000.0, 7000.0, 15000.0, 0.0, INFINITY, "Грозовой фронт"));
    
    // Запретная зона - цилиндрическая область
    optimizer->addConstraint(AirspaceConstraint::createCylinder(
        70000.0, 30000.0, 0.0, 9000.0, 10000.0, 0.0, INFINITY, "Запретная зона"));
    
    // Коридор для обязательного прохождения
    optimizer->addConstraint(AirspaceConstraint::createCorridor(
        30000.0, 20000.0, 8000.0, 60000.0, 50000.0, 9000.0, 5000.0, 0.0, INFINITY, "Обязательный коридор"));
    
    // Временное окно для прибытия
    optimizer->addConstraint(AirspaceConstraint::createTimeWindow(
        100000.0, 80000.0, 5000.0, 3600.0, 4500.0, "Слот прибытия"));
    
    // Оптимизация по расходу топлива
    std::cout << "\nОптимизация по расходу топлива:" << std::endl;
    TrajectoryResult result_fuel = optimizer->optimize4D("fuel", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 100000.0}, {"yf", 80000.0}, {"zf", 5000.0}
    });
    optimizer->visualizeTrajectory(result_fuel);
    optimizer->exportTrajectory(result_fuel, "fuel_trajectory.csv");
    
    // Оптимизация по времени
    std::cout << "\nОптимизация по времени:" << std::endl;
    TrajectoryResult result_time = optimizer->optimize4D("time", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 100000.0}, {"yf", 80000.0}, {"zf", 5000.0}
    });
    optimizer->visualizeTrajectory(result_time);
    optimizer->exportTrajectory(result_time, "time_trajectory.csv");
    
    // Сбалансированная оптимизация
    std::cout << "\nСбалансированная оптимизация (время + топливо):" << std::endl;
    TrajectoryResult result_balanced = optimizer->optimize4D("balanced", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 100000.0}, {"yf", 80000.0}, {"zf", 5000.0}
    });
    optimizer->visualizeTrajectory(result_balanced);
    optimizer->exportTrajectory(result_balanced, "balanced_trajectory.csv");
    
    // Сравнение результатов
    std::cout << "\n=== СРАВНЕНИЕ РЕЗУЛЬТАТОВ ===" << std::endl;
    std::cout << "Критерий   | Время полета (мин) | Расход топлива | Время вычислений (сек)" << std::endl;
    std::cout << "-----------|-------------------|---------------|----------------------" << std::endl;
    std::cout << "Топливо    | " << std::fixed << std::setprecision(2) 
              << result_fuel.total_flight_time / 60.0 << " | " 
              << result_fuel.fuel_consumption << " | " 
              << result_fuel.computation_time << std::endl;
    std::cout << "Время      | " << std::fixed << std::setprecision(2) 
              << result_time.total_flight_time / 60.0 << " | " 
              << result_time.fuel_consumption << " | " 
              << result_time.computation_time << std::endl;
    std::cout << "Баланс     | " << std::fixed << std::setprecision(2) 
              << result_balanced.total_flight_time / 60.0 << " | " 
              << result_balanced.fuel_consumption << " | " 
              << result_balanced.computation_time << std::endl;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    demoAPI();
    demo4DTrajectory();

    return 0;
}