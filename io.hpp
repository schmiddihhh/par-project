#include <mpi.h>
#include <cassert>
#include <string>
#include <regex>
#include <iostream>
#include <bitset>

void read_header( std::string filename, int* width_ptr, int* height_ptr, size_t* header_size_ptr, int* max_value_ptr );
bool* read_partial_image( std::string filename, size_t header_size, int width, int height, size_t offset, size_t stride, size_t pixel_size );
void write_partial_image( std::string filename, bool* partial_image, size_t header_size, int partition_width, int partition_height, size_t offset, size_t stride );
void write_image_header( std::string filename, int width, int height, size_t* header_size );
bool* read_file( std::string filename, int* width, int* height );
void write_file( std::string filename, bool* image, int width, int height );