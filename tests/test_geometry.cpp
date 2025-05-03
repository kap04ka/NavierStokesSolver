// #include "../src/core/Geometry.hpp"
// #include <cassert>
// #include <iostream>

// int main() {
//     constexpr std::size_t NX = 10, NY = 5;
//     cfd::Geometry geom(NX, NY, 1.0, 0.5);

//     // add small block in the centre
//     geom.add_rectangle(4,1,5,3);

//     geom.export_tags_csv("tags.csv");

//     // quick sanity check
//     assert(geom.tags()(4,2) == cfd::CellTag::FLUID);

//     std::cout << "Geometry smoke‑test passed.\n";
//     return 0;
// }