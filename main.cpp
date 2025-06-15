#define NOMINMAX

#include <iostream>
#include <windows.h>
#include <memory>
#include <chrono>
#include <vector>
#include <iomanip>
#include "TrajectoryOptimizer.hpp"

struct ScenarioResult {
    std::string name;
    std::string objective;
    std::string aircraft_type;
    int flight_time_min;
    int fuel_consumption_kg;
    int computation_time_ms;
    double distance_km;
    double avg_speed_kmh;
    bool has_constraints;
    bool success;
};

void runScenario(const std::string& name, std::unique_ptr<TrajectoryOptimizer>& optimizer,
                 const std::string& objective, const std::map<std::string, double>& target,
                 const std::string& aircraft_type, std::vector<ScenarioResult>& results) {
    
    std::cout << "🛩️  " << name << " ";
    
    auto start = std::chrono::high_resolution_clock::now();
    
    ScenarioResult scenario;
    scenario.name = name;
    scenario.objective = objective;
    scenario.aircraft_type = aircraft_type;
    scenario.has_constraints = !optimizer->getConstraints().empty();
    scenario.success = false;
    
    try {
        TrajectoryResult result = optimizer->optimize4D(objective, target);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        if (result.computation_time > 0) {
            scenario.flight_time_min = (int)(result.total_flight_time / 60.0);
            scenario.fuel_consumption_kg = (int)result.fuel_consumption;
            scenario.computation_time_ms = duration.count();
            scenario.success = true;
            
            // Вычисляем расстояние
            if (!result.positions.empty() && result.positions.size() >= 2) {
                auto& start_pos = result.positions[0];
                auto& end_pos = result.positions.back();
                scenario.distance_km = std::sqrt(
                    std::pow(end_pos[0] - start_pos[0], 2) + 
                    std::pow(end_pos[1] - start_pos[1], 2) + 
                    std::pow(end_pos[2] - start_pos[2], 2)
                ) / 1000.0;
                
                // Средняя скорость
                if (scenario.flight_time_min > 0) {
                    scenario.avg_speed_kmh = scenario.distance_km / (scenario.flight_time_min / 60.0);
                }
            }
            
            std::cout << "✅ Время: " << scenario.flight_time_min << " мин, "
                      << "Расход: " << scenario.fuel_consumption_kg << " кг, "
                      << "Вычисления: " << scenario.computation_time_ms << " мс" << std::endl;
            
            optimizer->exportTrajectory(result, name + "_trajectory.csv");
        } else {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            scenario.computation_time_ms = duration.count();
            std::cout << "❌ Неудача (" << duration.count() << " мс)" << std::endl;
        }
        
    } catch (std::exception& e) {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        scenario.computation_time_ms = duration.count();
        std::cout << "❌ Ошибка (" << duration.count() << " мс)" << std::endl;
    }
    
    results.push_back(scenario);
}

void printComparisonTable(const std::vector<ScenarioResult>& results) {
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "📊 СРАВНИТЕЛЬНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ" << std::endl;
    std::cout << std::string(120, '=') << std::endl;
    
    // Заголовок таблицы
    std::cout << std::left 
              << std::setw(20) << "Сценарий"
              << std::setw(12) << "Самолет"
              << std::setw(10) << "Цель"
              << std::setw(8) << "Время"
              << std::setw(10) << "Расход"
              << std::setw(10) << "Дистанция"
              << std::setw(12) << "Ср.скорость"
              << std::setw(12) << "Вычисления"
              << std::setw(12) << "Ограничения"
              << std::setw(8) << "Статус"
              << std::endl;
    
    std::cout << std::left 
              << std::setw(20) << ""
              << std::setw(12) << ""
              << std::setw(10) << ""
              << std::setw(8) << "(мин)"
              << std::setw(10) << "(кг)"
              << std::setw(10) << "(км)"
              << std::setw(12) << "(км/ч)"
              << std::setw(12) << "(мс)"
              << std::setw(12) << ""
              << std::setw(8) << ""
              << std::endl;
    
    std::cout << std::string(120, '-') << std::endl;
    
    // Данные таблицы
    for (const auto& result : results) {
        std::cout << std::left 
                  << std::setw(20) << result.name.substr(0, 19)
                  << std::setw(12) << result.aircraft_type
                  << std::setw(10) << result.objective;
        
        if (result.success) {
            std::cout << std::setw(8) << result.flight_time_min
                      << std::setw(10) << result.fuel_consumption_kg
                      << std::setw(10) << std::fixed << std::setprecision(0) << result.distance_km
                      << std::setw(12) << std::fixed << std::setprecision(0) << result.avg_speed_kmh
                      << std::setw(12) << result.computation_time_ms
                      << std::setw(12) << (result.has_constraints ? "Да" : "Нет")
                      << std::setw(8) << "✅";
        } else {
            std::cout << std::setw(8) << "-"
                      << std::setw(10) << "-"
                      << std::setw(10) << "-"
                      << std::setw(12) << "-"
                      << std::setw(12) << result.computation_time_ms
                      << std::setw(12) << (result.has_constraints ? "Да" : "Нет")
                      << std::setw(8) << "❌";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::string(120, '-') << std::endl;
    
    // Статистика
    int successful_runs = 0;
    int total_time = 0;
    int total_fuel = 0;
    double total_distance = 0;
    int total_computation = 0;
    
    for (const auto& result : results) {
        if (result.success) {
            successful_runs++;
            total_time += result.flight_time_min;
            total_fuel += result.fuel_consumption_kg;
            total_distance += result.distance_km;
        }
        total_computation += result.computation_time_ms;
    }
    
    std::cout << "\n📈 ОБЩАЯ СТАТИСТИКА:" << std::endl;
    std::cout << "• Успешных вычислений: " << successful_runs << "/" << results.size() 
              << " (" << (successful_runs * 100 / results.size()) << "%)" << std::endl;
    
    if (successful_runs > 0) {
        std::cout << "• Среднее время полета: " << (total_time / successful_runs) << " мин" << std::endl;
        std::cout << "• Средний расход топлива: " << (total_fuel / successful_runs) << " кг" << std::endl;
        std::cout << "• Общая дистанция: " << std::fixed << std::setprecision(0) << total_distance << " км" << std::endl;
    }
    
    std::cout << "• Общее время вычислений: " << (total_computation / 1000.0) << " сек" << std::endl;
    std::cout << "• Среднее время вычислений: " << (total_computation / results.size()) << " мс" << std::endl;
    
    // Найдем лучшие результаты
    if (successful_runs > 0) {
        auto fastest_flight = *std::min_element(results.begin(), results.end(),
            [](const ScenarioResult& a, const ScenarioResult& b) {
                return a.success && (!b.success || a.flight_time_min < b.flight_time_min);
            });
        
        auto most_efficient = *std::min_element(results.begin(), results.end(),
            [](const ScenarioResult& a, const ScenarioResult& b) {
                return a.success && (!b.success || a.fuel_consumption_kg < b.fuel_consumption_kg);
            });
        
        auto fastest_computation = *std::min_element(results.begin(), results.end(),
            [](const ScenarioResult& a, const ScenarioResult& b) {
                return a.success && (!b.success || a.computation_time_ms < b.computation_time_ms);
            });
        
        std::cout << "\n🏆 ЛУЧШИЕ РЕЗУЛЬТАТЫ:" << std::endl;
        std::cout << "• Самый быстрый полет: " << fastest_flight.name 
                  << " (" << fastest_flight.flight_time_min << " мин)" << std::endl;
        std::cout << "• Самый экономичный: " << most_efficient.name 
                  << " (" << most_efficient.fuel_consumption_kg << " кг)" << std::endl;
        std::cout << "• Самые быстрые вычисления: " << fastest_computation.name 
                  << " (" << fastest_computation.computation_time_ms << " мс)" << std::endl;
    }
    
    std::cout << std::string(120, '=') << std::endl;
}

void demonstrateFlightScenarios() {
    std::cout << "🚀 ДЕМОНСТРАЦИЯ АВИАЦИОННЫХ СЦЕНАРИЕВ\n" << std::endl;
    
    std::vector<ScenarioResult> results;
    
    // Параметры современного авиалайнера
    SystemParams aircraft(180000.0, 400000.0, 800000.0, 280.0, 1000.0, 13000.0);
    OptimizationParams optimization(20, 0.0, "ipopt");
    
    auto optimizer = std::make_unique<TrajectoryOptimizer>(aircraft, optimization);
    
    // 1. Короткий рейс (эффективность топлива)
    optimizer->clearConstraints();
    runScenario("short_haul_fuel", optimizer, "fuel", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 150000.0}, {"yf", 100000.0}, {"zf", 8000.0}
    }, "Авиалайнер", results);
    
    // 2. Средний рейс (минимальное время)
    runScenario("medium_haul_time", optimizer, "time", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 300000.0}, {"yf", 200000.0}, {"zf", 11000.0}
    }, "Авиалайнер", results);
    
    // 3. Дальний рейс (сбалансированная оптимизация)
    runScenario("long_haul_balanced", optimizer, "balanced", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 500000.0}, {"yf", 350000.0}, {"zf", 12000.0}
    }, "Авиалайнер", results);
    
    // 4. Высокогорный полет (энергоэффективность)
    runScenario("high_altitude_energy", optimizer, "energy", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 5000.0}, 
        {"xf", 250000.0}, {"yf", 150000.0}, {"zf", 12500.0}
    }, "Авиалайнер", results);
    
    // 5. Полет с одним препятствием
    optimizer->clearConstraints();
    optimizer->addConstraint(AirspaceConstraint::createSphere(
        125000.0, 75000.0, 7000.0, 15000.0, 0.0, INFINITY, "Запретная зона"));
    
    runScenario("obstacle_avoidance", optimizer, "balanced", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 3000.0}, 
        {"xf", 250000.0}, {"yf", 150000.0}, {"zf", 9000.0}
    }, "Авиалайнер", results);
    
    // 6. Транссонический полет
    SystemParams fighter(25000.0, 600000.0, 1200000.0, 600.0, 500.0, 18000.0);
    auto fighter_optimizer = std::make_unique<TrajectoryOptimizer>(fighter, optimization);
    fighter_optimizer->clearConstraints();
    
    runScenario("supersonic_flight", fighter_optimizer, "time", {
        {"x0", 0.0}, {"y0", 0.0}, {"z0", 8000.0}, 
        {"xf", 400000.0}, {"yf", 300000.0}, {"zf", 15000.0}
    }, "Истребитель", results);
    
    // Печатаем сравнительную таблицу
    printComparisonTable(results);
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    demonstrateFlightScenarios();

    return 0;
}