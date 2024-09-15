#include <iostream>
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


typedef std::vector<double> grid1;
typedef std::tuple<double, double, double> coordinates;
typedef std::tuple<double, double> vector_coordinates;
typedef std::vector<coordinates> grid3;



grid1 mesh (double left_border, const double & right_border, const double & mesh_step);

grid3 function_nodes (grid1 & x_grid, grid1 & y_grid, double f(double & x, double & y));

double kenny (double & x, double & y);

template<typename... Tp>
void file_creation (const std::string & file_name, std::vector<Tp...> & data);

void plot (const std::string & name, const std::string & title);

int main() {
    grid1 x_grid = std::move(mesh(x_min, x_max, step));
    grid1 y_grid = std::move(mesh(y_min, y_max, step));
    grid3 McCormick = std::move(function_nodes(x_grid, y_grid, kenny));
    file_creation ("kenny_test", McCormick);
    plot("kenny_test", "McCormick function");
    return 0;
}




double kenny (double & x, double & y) {
    return std::sin(x + y) + std::pow(x - y, 2.0) - 1.5*x + 2.5*y + 1.0;
}


grid3 function_nodes (grid1 & x_grid, grid1 & y_grid, double f(double & x, double & y)) {
    grid3 z;
    for (double & x : x_grid)
        for (double & y : y_grid)
            z.emplace_back(std::make_tuple(x, y, f(x, y)));
    return z;
}


bool is_equal (double & a, double b) {
    return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


double diff (double & f_prev, double & f_next, const double & q_step) {
    return (f_next - f_prev) / (2.0 * q_step);
}

// It should be rewrited is template if it would be usable as recursion of use array instead
void find_neighbours (grid3 & nodes, int & node, std::tuple<double, double> & x_neighbours, std::tuple<double, double> y_neighbours) {
    int y_left, y_rigth, x_left, x_right;
    //std::get<0>
    for (int i = node; i >= 0; --i)
        if (!is_equal(std::get<1>(nodes[i]), std::get<1>(nodes[node])) && is_equal(std::get<0>(nodes[i]), std::get<0>(nodes[node]))) {
            y_left = i;
            break;
        }
    for (int i = node; i < nodes.size(); ++i)
        if (!is_equal(std::get<1>(nodes[i]), std::get<1>(nodes[node])) && is_equal(std::get<0>(nodes[i]), std::get<0>(nodes[node]))) {
            y_rigth = i;
            break;
        }
    y_neighbours = std::make_tuple(std::get<2>(nodes[y_left]), std::get<2>(nodes[y_rigth]));
    //std::get<1>
    for (int i = node; i >= 0; --i)
        if (is_equal(std::get<1>(nodes[i]), std::get<1>(nodes[node])) && !is_equal(std::get<0>(nodes[i]), std::get<0>(nodes[node]))) {
            x_left = i;
            break;
        }
    for (int i = node; i < nodes.size(); ++i)
        if (is_equal(std::get<1>(nodes[i]), std::get<1>(nodes[node])) && !is_equal(std::get<0>(nodes[i]), std::get<0>(nodes[node]))) {
            x_right = i;
            break;
        }
    x_neighbours = std::make_tuple(std::get<2>(nodes[x_left]), std::get<2>(nodes[x_right]));
}


vector_coordinates grad (grid3 & nodes, int & node) {
    std::tuple<double, double> x_neighbours, y_neighbours;
    find_neighbours(nodes, node, x_neighbours, y_neighbours);
    double dfdx = //create template
    double dfdy = //!!
    return std::make_tuple
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