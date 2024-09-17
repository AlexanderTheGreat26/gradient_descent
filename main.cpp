#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <sstream>
#include <tuple>


const double x_min = -1.5;
const double x_max = 4.0;
const double y_min = -3.0;
const double y_max = 4.0;
const double step = 0.01;
const double eps = std::numeric_limits<double>::epsilon();
const double l = 0.1;

typedef std::vector<double> grid1;
//typedef std::array<double, 3> coordinates;
//typedef std::array<double, 2> vector_coordinates;
//typedef std::vector<coordinates> grid3;

typedef std::tuple<double, double, double> coordinates;
typedef std::tuple<double, double> vector_coordinates;
typedef std::vector<coordinates> grid3;


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


grid1 mesh (double left_border, const double & right_border, const double & mesh_step);

grid3 function_nodes (grid1 & x_grid, grid1 & y_grid, double f(double & x, double & y));

double kenny (double & x, double & y);

template<typename... Tp>
void file_creation (const std::string & file_name, std::vector<Tp...> & data);

void plot (const std::string & name, const std::string & title);

grid3 gradient_descent (coordinates q_start, coordinates neighbour_point, double f(double & x, double & y));




int main() {
    grid1 x_grid = std::move(mesh(x_min, x_max, step));
    grid1 y_grid = std::move(mesh(y_min, y_max, step));
    grid3 McCormick = std::move(function_nodes(x_grid, y_grid, kenny));
//    file_creation ("kenny_test", McCormick);
//    plot("kenny_test", "McCormick function");

    std::uniform_int_distribution<> dis_node (1, McCormick.size()-2);
    int start_node = 0; // = dis_node(gen);

    grid3 traj = gradient_descent(McCormick[start_node], McCormick[start_node+750], kenny);
    file_creation ("kenny_traj", traj);


    return 0;
}




double kenny (double & x, double & y) {
    return std::sin(x + y) + std::pow(x - y, 2.0) - 1.5*x + 2.5*y + 1.0;
    //return std::pow(x, 2) /10.0 + std::pow(y, 2);
}


grid3 function_nodes (grid1 & x_grid, grid1 & y_grid, double f(double & x, double & y)) {
    grid3 z;
    for (double & x : x_grid)
        for (double & y : y_grid)
            z.emplace_back(std::make_tuple(x, y, f(x, y)));
    return z;
}


double diff (double & f_curr, double & f_prev, double dx) {
    if (std::fabs(dx) < std::numeric_limits<double>::epsilon())
        return 0;
    else
        return (f_curr - f_prev) / dx;
}


vector_coordinates grad (vector_coordinates & current_point, vector_coordinates & new_point, double f (double & x, double & y)) {
    double x_curr = std::get<0>(current_point);
    double y_curr = std::get<1>(current_point);
    double x_new = std::get<0>(new_point);
    double y_new = std::get<1>(new_point);
    double dx = x_new - x_curr;
    double dy = y_new - y_curr;
    double f_curr = f(x_curr, y_curr);
    double f_x_new = f(x_new, y_curr);
    double f_y_new = f(x_curr, y_new);
    return std::make_tuple(diff(f_x_new, f_curr, dx), diff(f_y_new, f_curr, dy));
}


template<size_t Is = 0, typename... Tp>
void new_point (std::tuple<Tp...> & point, std::tuple<Tp...> & function_gradient, const double & lambda) {
        std::get<Is>(point) -= lambda * std::get<Is>(function_gradient);
    if constexpr (Is + 1 != sizeof...(Tp))
        new_point<Is + 1>(point, function_gradient, lambda);
}


//template<typename T>
grid3 gradient_descent (coordinates q_start, coordinates neighbour_point, double f(double & x, double & y)) {
    // We need to points to have numerical gradient. Let's take some random point and it's neighbour for grad calc.
    double old_value, new_value = std::get<2>(q_start);
    vector_coordinates current_point = std::make_tuple(std::get<0>(q_start), std::get<1>(q_start));
    vector_coordinates next_point = std::make_tuple(std::get<0>(neighbour_point), std::get<1>(neighbour_point));
    vector_coordinates grad_value = grad(current_point, next_point, f);

    grid3 trajectory {q_start, neighbour_point};



    do {
        int n = trajectory.size()-1; // length of trajectory.

        old_value = new_value;

        current_point = std::make_tuple(std::get<0>(trajectory[n]), std::get<1>(trajectory[n]));

        next_point = current_point;

        new_point(next_point, grad_value, l);

        grad_value = grad(current_point, next_point, f);

        new_value = f(std::get<0>(next_point), std::get<1>(next_point));

        trajectory.emplace_back(std::make_tuple(std::get<0>(next_point), std::get<1>(next_point), new_value));

    } while (std::fabs(old_value - new_value) > eps);
    return trajectory;
}





grid1 mesh (double left_border, const double & right_border, const double & mesh_step) {
    std::vector <double> xx ((right_border-left_border) / mesh_step + 1);
    xx[0] = left_border;
    std::generate(xx.begin()+1, xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


template<typename T, size_t... Is>
std::string tuple2string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple2string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple2string_impl(t, std::make_index_sequence<size>{});
}


template<typename... Tp>
void file_creation (const std::string & file_name, std::vector<Tp...> & data) {
    std::ofstream fout;
    fout.open(file_name, std::ios::trunc);
    for (auto & i : data)
        fout << tuple2string(i) << '\n';
    fout.close();
}



//template<typename... Tp>
//void file_creation (const std::string & file_name, grid3 & data) {
//    std::ofstream fout;
//    fout.open(file_name, std::ios::trunc);
//    int N = data[0].size();
//    for (auto & i : data) {
//        for (int j = 0; j < N; ++j)
//            fout << i[j] << '\t';
//        fout << '\n';
//    }
//    fout.close();
//}


void plot (const std::string & name, const std::string & title) {
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 1920, 1080 font \"Helvetica,30\"",
                                      "set output \'" + title + ".jpg\'",
                                      "set title \'" + title + R"(' font 'Helvetica,20')",
                                      "set grid xtics ytics ztics",
                                      "set autoscale xfix",
                                      "set key off",
                                      "set ticslevel 0",
                                      //"set border 4095",
                                      "set hidden3d",
                                      "set dgrid3d 50,50 qnorm 4",
                                      "splot \'" + name + ".txt\' w l",
                                      "set tics font \'Helvetica,14\'",
                                      R"(set xlabel 'x' font 'Helvetica,20')",
                                      R"(set ylabel 'y' font 'Helvetica,20')",
                                      R"(set zlabel 'z' font 'Helvetica,20')",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

//int period (const double & q_min, const double & q_max, const double & mesh_step) {
//    return (q_max - q_min) / mesh_step;
//}


//std::array<std::array<double,2>, 2> neighbours (grid3 & nodes, int & node, const double & mesh_step) {
//    int x_period = period(nodes[0][0], nodes[nodes.size()-1][0], mesh_step);
//    int y_period = period(nodes[0][1], nodes[nodes.size()-1][1], mesh_step);
//    std::array<double,2> x_neighbours = {nodes[2][node - y_period], nodes[node + y_period][2]};
//    std::array<double,2> y_neighbours = {nodes[2][node - x_period], nodes[node + x_period][2]};
//    return {x_neighbours, y_neighbours};
//}
//
//
//vector_coordinates grad (grid3 & nodes, int & node, const double & mesh_step) {
//    auto q_neighbours = neighbours(nodes, node, mesh_step);
//    double f_x_prev = q_neighbours[0][0];
//    double f_x_next = q_neighbours[0][1];
//    double f_y_prev = q_neighbours[1][0];
//    double f_y_next = q_neighbours[1][1];
//    double dfdx = diff(f_x_prev, f_x_next, mesh_step);
//    double dfdy = diff(f_y_prev, f_y_next, mesh_step);
//    return {dfdx, dfdy};
//}
