#ifndef EXPORT_UTILS_HPP
#define EXPORT_UTILS_HPP

#include "core/Field2D.hpp" // Убедись, что путь правильный
#include <fstream>
#include <string>
#include <vector>
#include <iostream> // Для cerr
#include <iomanip>  // Для setprecision

namespace cfd {

template<typename T>
void exportFieldToCSV(
    const Field2D<T>& field, 
    const std::string& filename,
    std::size_t i_start = 0, 
    std::size_t j_start = 0,
    std::size_t i_end = 0, // 0 означает до конца по умолчанию
    std::size_t j_end = 0   // 0 означает до конца по умолчанию
) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    std::size_t actual_i_end = (i_end == 0 || i_end >= field.nx()) ? field.nx() - 1 : i_end;
    std::size_t actual_j_end = (j_end == 0 || j_end >= field.ny()) ? field.ny() - 1 : j_end;
    
    std::size_t actual_i_start = std::max(static_cast<std::size_t>(0), i_start);
    std::size_t actual_j_start = std::max(static_cast<std::size_t>(0), j_start);

    if (actual_i_start > actual_i_end || actual_j_start > actual_j_end) {
        std::cerr << "Warning: Invalid region for export in " << filename 
                  << ". Region: i[" << actual_i_start << "-" << actual_i_end 
                  << "], j[" << actual_j_start << "-" << actual_j_end 
                  << "]. Field dims: " << field.nx() << "x" << field.ny() << std::endl;
        outfile.close();
        return;
    }
    
    outfile << std::fixed << std::setprecision(6);

    for (std::size_t j = actual_j_start; j <= actual_j_end; ++j) {
        for (std::size_t i = actual_i_start; i <= actual_i_end; ++i) {
            outfile << field(i, j);
            if (i < actual_i_end) {
                outfile << ",";
            }
        }
        outfile << "\n";
    }
    outfile.close();
}

} // namespace cfd

#endif // EXPORT_UTILS_HPP