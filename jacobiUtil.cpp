#include <jacobi.hpp>

void print_vector(const std::vector<int>& vec) {
    std::cout << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]";
    std::cout << std::endl;
}

void output_dat(const std::string& filename, DataPoint dp) {

    std::ofstream f;
    f.open(filename + ".csv", std::ios::app);

    if (!f.is_open()) {
        std::cerr << "Error: Could not open file " << filename + ".csv" << std::endl;
        return;
    }

    // Check if file is empty to write header
    f.seekp(0, std::ios::end);
    if (f.tellp() == 0) {
        // Write the header
        f << "n_cores,threads,n_grid_points,time_elapsed,error\n";
    }

    f << dp.n_cores << "," << dp.threads << "," << dp.n_grid_points << ","
      << dp.time_elapsed << "," << dp.error << "\n";
    
    f.close();

}