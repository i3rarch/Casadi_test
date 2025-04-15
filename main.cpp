#define NOMINMAX

#include <iostream>
#include <windows.h>  // Для SetConsoleOutputCP
#include <memory>
#include "TrajectoryOptimizer.hpp"

void testConstraint(const std::string& constraint_name, 
                   std::unique_ptr<TrajectoryOptimizer>& optimizer,
                   const AirspaceConstraint& constraint) {
    std::cout << "\nТестирование ограничения: " << constraint_name << std::endl;
    
    // Очистка предыдущих ограничений
    optimizer->clearConstraints();
    
    // Добавление тестового ограничения
    optimizer->addConstraint(constraint);
    
    // Сначала пробуем отладочную версию (упрощенную)
    std::cout << "Запуск отладочной оптимизации..." << std::endl;
    TrajectoryResult result_debug = optimizer->debug4D("fuel", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 100000.0}, {"yf", 80000.0}, {"zf", 5000.0}
    });
    
    if (result_debug.computation_time > 0) {
        optimizer->visualizeTrajectory(result_debug);
        optimizer->exportTrajectory(result_debug, "debug_" + constraint_name + ".csv");
    } else {
        std::cout << "Отладочная оптимизация не удалась." << std::endl;
    }
    
    // Пробуем полную оптимизацию с этим ограничением
    std::cout << "Запуск полной оптимизации..." << std::endl;
    try {
        TrajectoryResult result_full = optimizer->optimize4D("fuel", {
            {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
            {"xf", 100000.0}, {"yf", 80000.0}, {"zf", 5000.0}
        });
        
        optimizer->visualizeTrajectory(result_full);
        optimizer->exportTrajectory(result_full, "full_" + constraint_name + ".csv");
        
        // Сравнение результатов
        if (result_debug.computation_time > 0) {
            std::cout << "\nОграничение: " << constraint_name << " - СРАВНЕНИЕ" << std::endl;
            std::cout << "Метод     | Время (мин) | Топливо | Вычисление (сек)" << std::endl;
            std::cout << "----------|------------|---------|------------------" << std::endl;
            std::cout << "Отладка   | " << std::fixed << std::setprecision(2) 
                      << result_debug.total_flight_time / 60.0 << " | " 
                      << result_debug.fuel_consumption << " | " 
                      << result_debug.computation_time << std::endl;
            std::cout << "Полный    | " << std::fixed << std::setprecision(2) 
                      << result_full.total_flight_time / 60.0 << " | " 
                      << result_full.fuel_consumption << " | " 
                      << result_full.computation_time << std::endl;
        }
    } catch (std::exception& e) {
        std::cerr << "Оптимизация не удалась: " << e.what() << std::endl;
        std::cerr << "Попробуйте настроить параметры ограничений или параметры решателя." << std::endl;
    }
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
    
    // Тестирование сферического ограничения (положение изменено во избежание проблем с NaN)
    testConstraint("sphere", optimizer, 
        AirspaceConstraint::createSphere(65000.0, 25000.0, 7000.0, 15000.0, 0.0, INFINITY, "Грозовой фронт"));
    
    // Тестирование цилиндрического ограничения - перемещено для предотвращения прямого перелета
    // Радиус увеличен для облегчения обнаружения
    testConstraint("cylinder", optimizer,
        AirspaceConstraint::createCylinder(45000.0, 35000.0, 6000.0, 8000.0, 18000.0, 0.0, INFINITY, "Зона ограничения"));

    // Тестирование комбинированных ограничений - скорректированные позиции
    std::cout << "\nТестирование комбинированных ограничений..." << std::endl;
    optimizer->clearConstraints();
    optimizer->addConstraint(AirspaceConstraint::createSphere(65000.0, 25000.0, 7000.0, 15000.0, 0.0, INFINITY, "Грозовой фронт"));
    optimizer->addConstraint(AirspaceConstraint::createCylinder(45000.0, 35000.0, 6000.0, 8000.0, 18000.0, 0.0, INFINITY, "Зона ограничения"));
    
    try {
        TrajectoryResult result_combined = optimizer->optimize4D("fuel", {
            {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
            {"xf", 100000.0}, {"yf", 80000.0}, {"zf", 5000.0}
        });
        optimizer->visualizeTrajectory(result_combined);
        optimizer->exportTrajectory(result_combined, "combined_constraints.csv");
    } catch (std::exception& e) {
        std::cerr << "Комбинированная оптимизация не удалась: " << e.what() << std::endl;
    }
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    demo4DTrajectory();

    return 0;
}