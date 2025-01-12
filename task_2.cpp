#include "io.hpp"
#include <array>
#include <iostream>
#include <string>
#include <cassert>
#include <math.h>
#include <chrono>


void print_image( bool* image, std::string title, int width, int height )
{
    std::cout << title << std::endl;
    for( size_t y = 0; y < height; y++ ) {
        for( size_t x = 0; x < width; x++ ) {
            if( image[y * width + x] ) std::cout << "#";
            else std::cout << ".";
        }
        std::cout << std::endl;
    }
}


void parallel_game_of_life( std::string input_path, std::string output_path, int steps )
{
    // get the MPI rank information
    int comm_size, rank;
    MPI_Comm_size( MPI_COMM_WORLD,
                   &comm_size );
    MPI_Comm_rank( MPI_COMM_WORLD,
                   &rank ); 
    assert( (comm_size & (comm_size - 1)) == 0 );     // checks if comm_size is a power of 2

    // read the image metadata from the image file
    // we need these (at least the image width and height) to calculate a good partition of the image
    int image_width, image_height, max_value;
    size_t header_size;
    read_header( input_path, &image_width, &image_height, &header_size, &max_value );

    // we have to find a useful partition of the image into subgrids
    // we only allow powers of 2 for the number of processes since that makes partitioning much easier
    // the following loop calculates a partition where all parts have equal size (at least almost -> rounding)
    int x_partitions = 1, y_partitions = 1;
    for( int splits = 0; splits < log2(comm_size); splits++ )
    {
        if( image_width / x_partitions > image_height / y_partitions ) x_partitions++;
        else y_partitions++;
    } 

    // next, every process needs to find its own partition
    int x_partition_num = rank % x_partitions;
    int y_partition_num = rank / x_partitions;

    // calculate the width and height of the partitions
    float avg_partition_width = image_width / x_partitions;
    float avg_partition_height = image_height / y_partitions;

    // map the partition to coordinates in the complete image
    int x_start = x_partition_num * avg_partition_width;
    int y_start = y_partition_num * avg_partition_height;
    int x_end = x_partition_num == x_partitions - 1 ? image_width : (x_partition_num + 1) * avg_partition_width;
    int y_end = y_partition_num == y_partitions - 1 ? image_height : (y_partition_num + 1) * avg_partition_height;
    int partition_width = x_end - x_start;      // the value at x_end (or y_end, respectively) is excluded
    int partition_height = y_end - y_start;

    // find out which processes will calculate the adjacent partitions
    // we need this for communicating with them
    int top_rank = ( comm_size + rank - x_partitions ) % comm_size;
    int bottom_rank = ( rank + x_partitions ) % comm_size;
    int right_rank = y_partition_num * x_partitions + ( x_partition_num + 1 ) % x_partitions;
    int left_rank = y_partition_num * x_partitions + ( x_partition_num + (x_partitions - 1) ) % x_partitions;

    // read the part of the image that we need
    bool* partial_image;
    size_t pixel_size = max_value < 256 ? 1 : 2;
    partial_image = read_partial_image( input_path,
                                        header_size, 
                                        partition_width, 
                                        partition_height, 
                                        y_start * image_width + x_start, 
                                        image_width,
                                        pixel_size );
    
    // // debug: print the read image
    // if( rank == 0 ) print_image( partial_image, "before start image", partition_width, partition_height );

    // define the size of our partial image including the ghost cells
    int width = partition_width + 2, height = partition_height + 2;

    // declaring two arrays, so we can alternate between reading from one and writing to the other
    bool** images = new bool*[2];
    images[0] = new bool[width * height];
    images[1] = new bool[width * height];
    bool* read_from, * write_to;

    // copy the input image to one of our calculating arrays: images[0]
    for( size_t y = 0; y < partition_height; y++ )
    {
        for( size_t x = 0; x < partition_width; x++ )
        {
            images[0][(y + 1) * width + (x + 1)] = partial_image[y * partition_width + x];
        }
    }

    // // debug: print the read image
    // if( rank == 0 ) print_image( images[0], "start image", width, height );

    // these are some definitions we need for exchanging data with other ranks
    MPI_Request* top_bottom_send_requests = new MPI_Request[2];
    MPI_Request* top_bottom_recv_requests = new MPI_Request[2];
    MPI_Request* left_right_send_requests = new MPI_Request[2];
    MPI_Request* left_right_recv_requests = new MPI_Request[2];

    // we use an MPI vector to make sending columns easier
    MPI_Datatype column;
    MPI_Type_vector( height,
                     1,
                     width,
                     MPI_CXX_BOOL,
                     &column );
    MPI_Type_commit( &column );
    
    // now we will start the actual calculation time steps
    for( int i = 0; i < steps; i++ )
    {
        // find out which image is used for reading and which for writing in this iteration
        read_from = images[i % 2];
        write_to = images[(i + 1) % 2];

        // now do the communication with the surrounding ranks

        // first, we have to send the top and bottom image lines to the ranks above and underneath
        int status = MPI_Isend( &read_from[width + 1],  // send the first image line to the process above
                                width - 2,
                                MPI_CXX_BOOL,
                                top_rank,
                                i* 2 + 1,
                                MPI_COMM_WORLD,
                                &top_bottom_send_requests[0] );
        assert( status == MPI_SUCCESS );

        status = MPI_Isend( &read_from[(height - 2) * width + 1],   // send the last image line to the process underneath
                            width - 2,
                            MPI_CXX_BOOL,
                            bottom_rank,
                            i * 2,
                            MPI_COMM_WORLD,
                            &top_bottom_send_requests[1] );
        assert( status == MPI_SUCCESS );

        // now we will receive the lines from the ranks above and underneath
        status = MPI_Irecv( &read_from[1],  // receive the last image line from the process above and store it in the top ghost cells
                            width - 2,
                            MPI_CXX_BOOL,
                            top_rank,
                            i * 2,
                            MPI_COMM_WORLD,
                            &top_bottom_recv_requests[0] );
        assert( status == MPI_SUCCESS );

        status = MPI_Irecv( &read_from[ (height - 1) * width + 1],  // receive the first image line from the process underneath and store it in the bottom ghost cells
                            width - 2,
                            MPI_CXX_BOOL,
                            bottom_rank,
                            i * 2 + 1,
                            MPI_COMM_WORLD,
                            &top_bottom_recv_requests[1] );
        assert( status == MPI_SUCCESS );

        // wait for all reveices to complete
        status = MPI_Waitall( 2, top_bottom_recv_requests, MPI_STATUSES_IGNORE );

        // now, all ghost cells directly above and underneath the image are filled
        // all ghost cells to the left and right of the image, including the corners, are still missing
        // these will be exchanged now

        status = MPI_Isend( &read_from[width - 2],  // send the rightmost image column (including corner ghost cells) to the right process
                            1,
                            column,
                            right_rank,
                            i* 2 + 1,
                            MPI_COMM_WORLD,
                            &left_right_send_requests[0] );
        assert( status == MPI_SUCCESS );

        status = MPI_Isend( &read_from[1],   // send the leftmost image column (incl. corner ghost cells) to the left process
                            1,
                            column,
                            left_rank,
                            i * 2,
                            MPI_COMM_WORLD,
                            &left_right_send_requests[1] );
        assert( status == MPI_SUCCESS );

        // now we will receive the lines from the ranks th the left and right
        status = MPI_Irecv( &read_from[width - 1],  // receive the leftmost image column (incl. corner ghost cells) from the right process
                            1,
                            column,
                            right_rank,
                            i * 2,
                            MPI_COMM_WORLD,
                            &left_right_recv_requests[0] );
        assert( status == MPI_SUCCESS );

        status = MPI_Irecv( &read_from[0],  // receive the rightmost image column (incl. corner ghost cells) from the left process
                            1,
                            column,
                            left_rank,
                            i * 2 + 1,
                            MPI_COMM_WORLD,
                            &left_right_recv_requests[1] );
        assert( status == MPI_SUCCESS );

        // wait for all receives to complete
        status = MPI_Waitall( 2, left_right_recv_requests, MPI_STATUSES_IGNORE );
        assert( status == MPI_SUCCESS );

        // now apply the stencil operation to all grid points
        for( int y = 1; y < height - 1; y++ )
        {
            for( int x = 1; x < width - 1; x++ )
            {
                int living_cells = 0;
                for( int ydiff = -1; ydiff <= 1; ydiff++ )
                {
                    living_cells += read_from[(y + ydiff) * width + (x - 1)];
                    living_cells += read_from[(y + ydiff) * width + x];
                    living_cells += read_from[(y + ydiff) * width + (x + 1)];
                }
                write_to[y * width + x] = living_cells == 3 || (read_from[y * width + x] == 1 && living_cells == 4);
            }
        }
    }

    // finally, we have to write the result to the given output path
    // rank 0 writes the image header and broadcasts its size to all other processes
    // (all processes need this size to calculate the offset to write at)
    if( rank == 0 )
    {
        write_image_header( output_path,
                            image_width,
                            image_height,
                            &header_size );
    }

    int status = MPI_Bcast( &header_size,
                            1,
                            MPI_INT,
                            0,
                            MPI_COMM_WORLD );
    assert( status == MPI_SUCCESS );

    bool* calc_result = images[steps % 2];
    for( int y = 0; y < partition_height; y++ )
    {
        for( int x = 0; x < partition_width; x++ )
        {
            partial_image[y * partition_width + x] = calc_result[(y + 1) * width + (x + 1)];
        }
    }
    write_partial_image( output_path,
                         partial_image,
                         header_size,
                         partition_width,
                         partition_height,
                         y_start * image_width + x_start,
                         image_width );

    // // erstmal debug print
    // if( rank == 0 ) print_image( final_image, "final image", width, height );
    
    delete[] top_bottom_send_requests;
    delete[] top_bottom_recv_requests;
    delete[] left_right_send_requests;
    delete[] left_right_recv_requests;
    delete[] images[0];
    delete[] images[1];
    delete[] images;
}


void run_test( std::string testname, std::string input_file, std::string output_file, std::string expected_output_file, int time_steps, bool debug_print = false )
{
    std::cout << "starting test: " << testname << std::endl;

    // read the test files
    bool* input, * output, * expected_outcome;
    int width, height;
    input = read_file(input_file, &width, &height);
    expected_outcome = read_file(expected_output_file, &width, &height);

    // debug: print the read files on the console
    if( debug_print ) 
    {
        print_image(input, "INPUT", width, height);
        print_image(expected_outcome, "EXPECTED", width, height);
    }

    // execute the game of life on the input and read the resulting image
    parallel_game_of_life( input_file, output_file, time_steps );
    output = read_file( output_file, &width, &height );

    // debug: print the resulting image to the console
    if( debug_print ) print_image(output, "RESULT", width, height);

    // compare the result to the expected outcome
    int test_passed = true;
    for( size_t i = 0; i < width * height; i++ )
    {
        if( ! output[i] == expected_outcome[i] )
        {
            std::cout << " !!> TEST FAILED: " << testname << std::endl;
            test_passed = false;
            break;
        }
    }
    if( test_passed ) std::cout << " --> test passed: " << testname << std::endl;

    delete[] input;
    delete[] output;
    delete[] expected_outcome;
}


void test_game_of_life( bool debug_print = false )
{
    std::cout << "---- testing the game_of_life function ----" << std::endl;

    // test 1: a cell should die if it has less than 2 or more than 3 neighbors
    run_test("dying cell", "test_images/dying_cell_input.ppm", "test_images/test_out.ppm", "test_images/dying_cell_1_step.ppm", 1, debug_print);

    // test 2: a cell should stay alive if it has 2 or 3 neighbors
    run_test("surviving cell", "test_images/surviving_cell_input.ppm", "test_images/test_out.ppm", "test_images/surviving_cell_1_step.ppm", 1, debug_print);

    // test 3: a cell should resurrect if it has exactly 3 neighbors
    run_test("resurrecting cell", "test_images/resurrecting_cell_input.ppm", "test_images/test_out.ppm", "test_images/resurrecting_cell_1_step.ppm", 1, debug_print);

    // test 4: glider
    run_test("glider", "test_images/glider_input.ppm", "test_images/test_out.ppm", "test_images/glider_4_steps.ppm", 4, debug_print);

    // test 5: periodic boundary (glider goes from bottom right corner to top left corner)
    run_test("periodic boundary", "test_images/periodic_boundary_input.ppm", "test_images/test_out.ppm", "test_images/periodic_boundary_n_steps.ppm", 24, debug_print);

    // test 6: spacefiller
    run_test("spacefiller", "test_images/spacefiller_input.ppm", "test_images/test_out.ppm", "test_images/spacefiller_150_steps.ppm", 150, debug_print);

    std::cout << "---- finished testing the game_of_life function ----" << std::endl;
}


int main()
{
    MPI_Init(NULL, NULL);

    test_game_of_life();

    // for( int i = 0; i < 100; i++ )
    // {
    //     std::string writename = "test_images/arschlecken" + std::to_string(i) + ".ppm";
    //     parallel_game_of_life( "test_images/spacefiller_input.ppm", writename, i);

    // }
    int iterations = 100;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    parallel_game_of_life( "test_images/spacefiller_input.ppm", "test_images/output.ppm", iterations);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "abs difference: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " milliseconds" << std::endl;
    std::cout << "time per iteration: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / iterations << "microseconds" << std::endl;

    MPI_Finalize();
    return 0;
}
