#include <iostream>
#include <casadi/casadi.hpp>
#include <windows.h>  // Add this include for SetConsoleOutputCP

using namespace casadi;

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    // Создаем символическую переменную (3-мерный вектор)
    MX x = MX::sym("x", 3);
    
    // Создаем трехмерную функцию (например, "чаша Розенброка" с дополнительным измерением)
    // f(x,y,z) = 100*(y-x^2)^2 + (1-x)^2 + (z-y^2)^2 + (1-y)^2
    MX f = 100*pow(x(1)-pow(x(0),2),2) + pow(1-x(0),2) + 
           100*pow(x(2)-pow(x(1),2),2) + pow(1-x(1),2);
    
    // Создаем объект для решения задачи оптимизации
    Dict opts;
    opts["ipopt.print_level"] = 0;  // Уменьшаем вывод информации от решателя
    opts["print_time"] = 0;
    MXDict nlp = {{"x", x}, {"f", f}};
    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    
    // Задаем начальную точку и границы
    std::vector<double> x0 = {-1.0, -1.0, -1.0};  // Начальная точка
    std::vector<double> lbx = {-10.0, -10.0, -10.0};  // Нижние границы
    std::vector<double> ubx = {10.0, 10.0, 10.0};     // Верхние границы
    
    // Формируем входные аргументы для решателя
    DMDict arg = {
        {"x0", DM(x0)},
        {"lbx", DM(lbx)},
        {"ubx", DM(ubx)}
    };
    
    // Решаем задачу оптимизации
    DMDict res = solver(arg);
    
    // Получаем оптимальное решение
    DM x_opt = res.at("x");
    DM f_opt = res.at("f");
    
    // Выводим результат
    std::cout << "Оптимальное решение x* = " << x_opt << std::endl;
    std::cout << "Минимальное значение f(x*) = " << f_opt << std::endl;
    
    return 0;
}