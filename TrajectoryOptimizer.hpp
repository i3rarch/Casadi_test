#pragma once

#include <casadi/casadi.hpp>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>

using namespace casadi;

class TrajectoryOptimizer {
private:
    // Параметры системы
    double mass;     // масса объекта
    double dt;       // шаг по времени
    int N;           // количество шагов горизонта
    
public:
    TrajectoryOptimizer(double mass = 1.0, double dt = 0.1, int N = 50);
    
    // Метод для построения оптимальной траектории
    DMDict trajectory(const std::string& objective = "fuel", 
                    const std::map<std::string, double>& kwargs = {});
    
    // Метод для визуализации результатов
    void visualizeTrajectory(const DMDict& result);
};