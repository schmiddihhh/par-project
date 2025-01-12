#include "io.hpp"


void read_header( std::string filename, int* width_ptr, int* height_ptr, size_t* header_size_ptr, int* max_value_ptr )
{
    // read the first characters of the file to get information about the file size
    MPI_File inputfile;
    int status = MPI_File_open( MPI_COMM_WORLD,
                                filename.c_str(),
                                MPI_MODE_RDONLY,
                                MPI_INFO_NULL,
                                &inputfile );
    assert( status == MPI_SUCCESS );
    char* fileheader = new char[100];
    status = MPI_File_read( inputfile,
                            fileheader,
                            100,
                            MPI_CHAR,
                            MPI_STATUS_IGNORE );
    assert( status == MPI_SUCCESS );
    status = MPI_File_close( &inputfile );
    assert( status == MPI_SUCCESS );

    // use a regex to extract the values: magic number, width and height of the image, maximum color value
    std::string s(fileheader);
    std::regex ppm_file(R"((?:#[^\n\r]*[\n\r]+)?(P6)\s+(?:#[^\n\r]*[\n\r]+)?(\d+)\s+(\d+)\s+(?:#[^\n\r]*[\n\r]+)?(\d+)\s+(?:#[^\n\r]*[\n\r]+)?)");
    std::smatch matches;
    std::regex_search(s, matches, ppm_file);
    
    *width_ptr = std::stoi(matches.str(2));
    *height_ptr = std::stoi(matches.str(3));
    *max_value_ptr = std::stoi(matches.str(4));
    *header_size_ptr = matches.str(0).length();

    delete[] fileheader;

    return;
}


bool* read_partial_image( std::string filename, size_t header_size, int partition_width, int partition_height, size_t offset, size_t stride, size_t pixel_size_per_color )
{
    /*
    Arguments:
    filename - path to the file that should be read
    header_size - the length of the PPM fileheader in chars
    partition_width - the width of the image
    partition_height - the height of the image
    offset - offset in pixels
    stride - stride in pixels
    pixel_size_per_color - size of the color values of each pixel in bytes
    */

    // map the logical size of the image to the actual file size
    size_t raw_image_size = 3 * pixel_size_per_color * partition_width * partition_height;     // number of colors per pixel (RGB) * size of each color value * number of pixels
    MPI_Offset raw_offset = header_size + 3 * pixel_size_per_color * offset;
    size_t raw_partition_width = 3 * pixel_size_per_color * partition_width;
    size_t raw_stride = 3 * pixel_size_per_color * stride;

    // create an array to hold the raw image data
    char* raw_image = new char[raw_image_size];

    // create an MPI datatype to read a partition from the whole image file
    MPI_Datatype partition;
    MPI_Type_vector( partition_height, raw_partition_width, raw_stride, MPI_CHAR, &partition );
    MPI_Type_commit( &partition );

    // read the partition into the array
    MPI_File inputfile;
    int status = MPI_File_open( MPI_COMM_WORLD,
                                filename.c_str(),
                                MPI_MODE_RDONLY,
                                MPI_INFO_NULL,
                                &inputfile );
    assert( status == MPI_SUCCESS );
    status = MPI_File_set_view( inputfile,
                                raw_offset,
                                MPI_CHAR,
                                partition,
                                "native",
                                MPI_INFO_NULL );
    assert( status == MPI_SUCCESS );
    status = MPI_File_read( inputfile,
                            raw_image,
                            raw_image_size,
                            MPI_CHAR,
                            MPI_STATUS_IGNORE );
    assert( status == MPI_SUCCESS );
    status = MPI_File_close( &inputfile );
    assert( status == MPI_SUCCESS );

    // read the image data from the file
    size_t pixel_count = partition_width * partition_height;
    size_t pixel_spacing = 3 * pixel_size_per_color;

    bool* image = new bool[pixel_count];
    char* position = raw_image;
    if( pixel_size_per_color == 1 )
    {
        for( int i = 0; i < pixel_count; i++ )
        {
            image[i] = *position == 0;
            position += pixel_spacing;
        }
    }
    else
    {
        for( int i = 0; i < pixel_count; i++ )
        {
            image[i] = (int16_t)(*position) << 8 + *(position + 1) == 0;
            position += pixel_spacing;
        }
    }

    return image;
}


void write_partial_image( std::string filename, bool* partial_image, size_t header_size, int partition_width, int partition_height, size_t offset, size_t stride )
{
    /*
    Arguments:
    filename - path to the file that should be read
    partial_image - bool array with the image values
    header_size - the length of the PPM fileheader in chars
    width - the width of the image
    height - the height of the image
    offset - offset in pixels
    stride - stride in pixels
    */

    // map the logical size of the image to the actual file size
    size_t partition_size = partition_width * partition_height;
    size_t raw_partition_size = 3 * partition_size;     // number of colors per pixel (RGB) * size of each color value * number of pixels
    size_t raw_partition_width = 3 * partition_width;
    size_t raw_stride = 3 * stride;
    size_t raw_offset = header_size + 3 * offset;

    // create an array to hold the raw image data
    char* raw_partition = new char[raw_partition_size];

    // convert the bool array to raw image data
    for( int i = 0; i < partition_size; i++ )
    {
        raw_partition[i * 3] = ! partial_image[i];
        raw_partition[i * 3 + 1] = ! partial_image[i];
        raw_partition[i * 3 + 2] = ! partial_image[i];
    }

    // create an MPI datatype to write the raw partition to the file
    MPI_Datatype partition;
    MPI_Type_vector( partition_height, raw_partition_width, raw_stride, MPI_CHAR, &partition );
    MPI_Type_commit( &partition );

    // write the 
    MPI_File outputfile;
    int status = MPI_File_open( MPI_COMM_WORLD,
                                filename.c_str(),
                                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                MPI_INFO_NULL,
                                &outputfile );
    assert( status == MPI_SUCCESS );
    status = MPI_File_set_view( outputfile,
                                raw_offset,
                                MPI_CHAR,
                                partition,
                                "native",
                                MPI_INFO_NULL );
    assert( status == MPI_SUCCESS );
    status = MPI_File_write( outputfile,
                             raw_partition,
                             raw_partition_size,
                             MPI_CHAR,
                             MPI_STATUS_IGNORE );
    assert( status == MPI_SUCCESS );
    status = MPI_File_close( &outputfile );
    assert( status == MPI_SUCCESS );

    return;
}


void write_image_header( std::string filename, int width, int height, size_t* header_size )
{
    // construct the ppm file header
    std::string ppm_header = "P6\n";    // magic number for binary ppm
    ppm_header += std::to_string(width) + " " + std::to_string(height) + '\n';  // width and height of the image
    ppm_header += "1\n";  // max pixel value

    *header_size = ppm_header.length();

    // write the header to the file
    MPI_File outputfile;
    int status = MPI_File_open( MPI_COMM_SELF,
                                filename.c_str(),
                                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                MPI_INFO_NULL,
                                &outputfile );
    assert( status == MPI_SUCCESS );
    status = MPI_File_write( outputfile,
                             ppm_header.c_str(),
                             ppm_header.length(),
                             MPI_CHAR,
                             MPI_STATUS_IGNORE );
    assert( status == MPI_SUCCESS );
    status = MPI_File_close( &outputfile );
    assert( status == MPI_SUCCESS );

    return;
}


bool* read_file( std::string filename, int* width_ptr, int* height_ptr )
{
    // open the file
    MPI_File inputfile;
    int status = MPI_File_open( MPI_COMM_WORLD,
                                filename.c_str(),
                                MPI_MODE_RDONLY,
                                MPI_INFO_NULL,
                                &inputfile );
    assert( status == MPI_SUCCESS );

    // determine the file size
    MPI_Offset filesize;
    status = MPI_File_get_size( inputfile,
                                &filesize );
    assert( status == MPI_SUCCESS );
    
    // read the file into an array
    char* input = new char[filesize];
    status = MPI_File_read_all( inputfile,
                                input,
                                filesize,
                                MPI_CHAR,
                                MPI_STATUS_IGNORE );
    assert( status == MPI_SUCCESS );

    // close the file
    status = MPI_File_close( &inputfile );
    assert( status == MPI_SUCCESS );

    // read the metadata from the file: magic number, width and height of the image, maximum color value
    std::string s(input);
    std::regex ppm_file(R"((?:#[^\n\r]*[\n\r]+)?(P6)\s+(?:#[^\n\r]*[\n\r]+)?(\d+)\s+(\d+)\s+(?:#[^\n\r]*[\n\r]+)?(\d+)\s+(?:#[^\n\r]*[\n\r]+)?)");
    std::smatch matches;
    std::regex_search(s, matches, ppm_file);
    
    std::string magic_number = matches.str(1);
    *width_ptr = std::stoi(matches.str(2));
    *height_ptr = std::stoi(matches.str(3));
    int maxval = std::stoi(matches.str(4));

    // read the image data from the file
    size_t pixel_size = maxval < 256 ? 1 : 2;
    size_t pixel_count = *width_ptr * *height_ptr;
    size_t pixel_spacing = 3 * pixel_size;

    bool* image = new bool[pixel_count];
    char* position = input + matches.str(0).length();
    if( pixel_size == 1 )
    {
        for( int i = 0; i < pixel_count; i++ )
        {
            image[i] = *position == 0;
            position += pixel_spacing;
        }
    }
    else
    {
        for( int i = 0; i < pixel_count; i++ )
        {
            image[i] = (int16_t)(*position) << 8 + *(position + 1) == 0;
            position += pixel_spacing;
        }
    }
    

    // free all arrays
    delete[] input;

    return image;
}


void write_file( std::string filename, bool* image, int width, int height )
{
    // create the ppm file header as a string
    std::string ppm_header = "P6\n";    // magic number for binary ppm
    ppm_header += std::to_string(width) + " " + std::to_string(height) + '\n';  // width and height of the image
    ppm_header += "1\n";  // max pixel value

    // create a char array to fit the whole ppm file
    int filesize = ppm_header.length() + width * height * 3;
    char* ppm_file = new char[filesize];

    // write the header to the char array
    for( int i = 0; i < ppm_header.length(); i++ )
    {
        ppm_file[i] = ppm_header[i];
    }

    // write the pixel values to the char array
    char* ppm_values = ppm_file + ppm_header.length();
    for( int i = 0; i < width * height; i++ )
    {
        ppm_values[i * 3] = ! image[i];
        ppm_values[i * 3 + 1] = ! image[i];
        ppm_values[i * 3 + 2] = ! image[i];
    }

    // open the file
    MPI_File outputfile;
    int status = MPI_File_open( MPI_COMM_WORLD,
                                filename.c_str(),
                                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                MPI_INFO_NULL,
                                &outputfile );
    assert( status == MPI_SUCCESS );

    // write all data to the file at once
    status = MPI_File_write( outputfile,
                             ppm_file,
                             filesize,
                             MPI_CHAR,
                             MPI_STATUS_IGNORE );
    assert( status == MPI_SUCCESS );

    // close the file
    status = MPI_File_close( &outputfile );
    assert( status == MPI_SUCCESS );

    // free all arrays
    delete[] ppm_file;

    return;
}