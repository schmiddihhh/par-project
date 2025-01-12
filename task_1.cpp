#include "io.hpp"
#include <array>
#include <iostream>
#include <string>
#include <cassert>
#include <chrono>


bool* game_of_life( int image_width, int image_height, bool *image, int steps, bool measure_update_time = false )
{
    /*
    The image array is saved row major.
    */

    // first, add a frame with a width of one pixel around the image to store the ghost cells
    int width = image_width + 2, height = image_height + 2;

    // declaring two arrays, so we can alternate between reading from one and writing to the other
    bool** images = new bool*[2];
    images[0] = new bool[width * height];
    images[1] = new bool[width * height];
    bool* read_from, * write_to;

    // copy the input image to one of our calculating arrays: images[0]
    for( size_t y = 0; y < image_height; y++ )
    {
        for( size_t x = 0; x < image_width; x++ )
        {
            images[0][(y + 1) * width + (x + 1)] = image[y * image_width + x];
        }
    }

    // measure the starting time for the benchmark
    std::chrono::steady_clock::time_point benchmark_begin, benchmark_end;
    if( measure_update_time )
    {
        benchmark_begin = std::chrono::steady_clock::now();
        std::cout << "\n---- Benchmark of game_of_life ----" << std::endl;
        std::cout << "starting time measurement" << std::endl;
    }
    for( int i = 0; i < steps; i++ )
    {
        // find out which image is used for reading and which for writing in this iteration
        read_from = images[i % 2];
        write_to = images[(i + 1) % 2];

        // first, copy the outer lines to the ghost cells
        for( int x = 1; x < width - 1; x++ )
        {
            // copy the top image line to the bottom ghost cells
            read_from[(height - 1) * width + x] = read_from[1 * width + x];
            // copy the bottom image line to the top ghost cells
            read_from[x] = read_from[(height - 2) * width + x];
        }
        // copy the outer columns to the ghost cells
        // include the corner ghost cells
        for( int y = 0; y < height; y++ )
        {
            // copy the left image column to the right ghost cells
            read_from[y * width + (width - 1)] = read_from[y * width + 1];
            // copy the right image column to the left ghost cells
            read_from[y * width] = read_from[y * width + (width - 2)];
        }
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

    // measure the finish time and calculate the avg time per iteration
    if( measure_update_time )
    {
        benchmark_end = std::chrono::steady_clock::now();
        std::cout << "finished time measurement" << std::endl;
        std::cout << "total elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(benchmark_end - benchmark_begin).count() << " milliseconds" << std::endl;
        std::cout << "time per iteration: " << std::chrono::duration_cast<std::chrono::microseconds>(benchmark_end - benchmark_begin).count() / steps << " microseconds" << std::endl;
    }

    // finally, we have to copy the result to a new array and return it
    bool* result = new bool[image_width * image_height];
    for( size_t y = 0; y < image_height; y++ )
    {
        for( size_t x = 0; x < image_width; x++ )
        {
            result[y * image_width + x] = images[steps % 2][(y + 1) * width + (x + 1)];
        }
    }

    delete[] images[0];
    delete[] images[1];
    delete[] images;

    return result;
}


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


void run_test( std::string testname, std::string input_file, std::string expected_output_file, int time_steps, bool debug_print = false )
{
    std::cout << "starting test: " << testname << std::endl;

    // read the test files
    bool* input, * expected_outcome;
    int width, height;
    input = read_file(input_file, &width, &height);
    expected_outcome = read_file(expected_output_file, &width, &height);

    // debug: print the read files on the console
    if( debug_print ) 
    {
        print_image(input, "INPUT", width, height);
        print_image(expected_outcome, "EXPECTED", width, height);
    }

    // execute the game of life on the input
    bool* result = game_of_life(width, height, input, time_steps);

    // debug: print the resulting image to the console
    if( debug_print ) print_image(result, "RESULT", width, height);

    // compare the result to the expected outcome
    int test_passed = true;
    for( size_t i = 0; i < width * height; i++ )
    {
        if( ! result[i] == expected_outcome[i] )
        {
            std::cout << " !!> TEST FAILED: " << testname << std::endl;
            test_passed = false;
            break;
        }
    }
    if( test_passed ) std::cout << " --> test passed: " << testname << std::endl;
}


void test_game_of_life( bool debug_print = false )
{
    std::cout << "---- testing the game_of_life function ----" << std::endl;

    // test 1: a cell should die if it has less than 2 or more than 3 neighbors
    run_test("dying cell", "test_images/dying_cell_input.ppm", "test_images/dying_cell_1_step.ppm", 1, debug_print);

    // test 2: a cell should stay alive if it has 2 or 3 neighbors
    run_test("surviving cell", "test_images/surviving_cell_input.ppm", "test_images/surviving_cell_1_step.ppm", 1, debug_print);

    // test 3: a cell should resurrect if it has exactly 3 neighbors
    run_test("resurrecting cell", "test_images/resurrecting_cell_input.ppm", "test_images/resurrecting_cell_1_step.ppm", 1, debug_print);

    // test 4: glider
    run_test("glider", "test_images/glider_input.ppm", "test_images/glider_4_steps.ppm", 4, debug_print);

    // test 5: periodic boundary (glider goes from bottom right corner to top left corner)
    run_test("periodic boundary", "test_images/periodic_boundary_input.ppm", "test_images/periodic_boundary_n_steps.ppm", 24, debug_print);

    // test 6: spacefiller
    run_test("spacefiller", "test_images/spacefiller_input.ppm", "test_images/spacefiller_150_steps.ppm", 150, debug_print);

    std::cout << "---- finished testing the game_of_life function ----" << std::endl;
}


int main()
{
    MPI_Init(NULL, NULL);
    
    // running all tests to make sure that game_of_life works correctly
    test_game_of_life();

    // benchmark the function by updating a 25'600 x 12'800 grid
    int benchmark_width = 25600, benchmark_height = 12800;
    bool* benchmark_input = new bool[benchmark_width * benchmark_height];
    game_of_life( benchmark_width, benchmark_height, benchmark_input, 10, true );

    MPI_Finalize();
    return 0;
}
