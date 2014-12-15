/*    
 * Copyright (c) 2014 K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <sstream>

#include <catch.hpp>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include "Atom/printFunctions.hpp"

namespace atom
{
namespace tests
{

int computeCartesianToTwoLineElementsTestFunction( 
    const gsl_vector* independentVariables, void* parameters, gsl_vector* residuals )
{
    gsl_vector_set( residuals, 0, 1.2 );
    gsl_vector_set( residuals, 1, 2.3 );
    gsl_vector_set( residuals, 2, 3.4 );
    gsl_vector_set( residuals, 3, 4.5 );
    gsl_vector_set( residuals, 4, 5.6 );
    gsl_vector_set( residuals, 5, 6.7 );
    return GSL_SUCCESS;
}

int computeAtomTestFunction( 
    const gsl_vector* independentVariables, void* parameters, gsl_vector* residuals )
{
    gsl_vector_set( residuals, 0, 1.2 );
    gsl_vector_set( residuals, 1, 2.3 );
    gsl_vector_set( residuals, 2, 3.4 );
    return GSL_SUCCESS;
}

struct Parameters{ };

TEST_CASE( "Print Cartesian state to Two-Line-Elements table header", "[print]" )
{
    // Set expected output string.    
    std::ostringstream tableHeader;
    tableHeader
        << "#  a              e              i              AoP            RAAN           "
        << "TA             f1             f2             f3             f4             "
        << "f5             f6             "
        << std::endl;

    REQUIRE( printCartesianToTleSolverStateTableHeader( ) == tableHeader.str( ) );
}

TEST_CASE( "Print Cartesian state to Two-Line-Elements solver state", "[print]" )
{
    // Set initial guess for GSL solver.
    gsl_vector* initialGuess = gsl_vector_alloc( 6 );
    gsl_vector_set( initialGuess, 0, 0.1 );
    gsl_vector_set( initialGuess, 1, 0.2 );      
    gsl_vector_set( initialGuess, 2, 0.3 );      
    gsl_vector_set( initialGuess, 3, 0.4 );      
    gsl_vector_set( initialGuess, 4, 0.5 );      
    gsl_vector_set( initialGuess, 5, 0.6 );      

    // Set up dummy GSL solver.
    Parameters parameters;

    gsl_multiroot_function testFunction 
        = { &computeCartesianToTwoLineElementsTestFunction, 
            6, 
            &parameters };

    const gsl_multiroot_fsolver_type* solverType = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc( solverType, 6 );
    gsl_multiroot_fsolver_set( solver, &testFunction, initialGuess );

    // Set expected output string.
    std::ostringstream tableRow;
    tableRow << "0  0.1            0.2            0.3            0.4            0.5            "
             << "0.6            1.2            2.3            3.4            4.5            "
             << "5.6            6.7            "
             << std::endl;

    REQUIRE( printCartesianToTleSolverState( 0, solver ) == tableRow.str( ) );
}

TEST_CASE( "Print Atom solver table header", "[print-atom-solver]")
{
    // Set expected output string.    
    std::ostringstream tableHeader;
    tableHeader
        << "#  v1_x           v1_y           v1_z           "
        << "f1             f2             f3             "
        << std::endl;

    REQUIRE( printAtomSolverStateTableHeader( ) == tableHeader.str( ) );
}

TEST_CASE( "Print Atom solver state", "[print]" )
{
    // Set initial guess for GSL solver.
    gsl_vector* initialGuess = gsl_vector_alloc( 3 );
    gsl_vector_set( initialGuess, 0, 0.1 );
    gsl_vector_set( initialGuess, 1, 0.2 );      
    gsl_vector_set( initialGuess, 2, 0.3 );     

    // Set up dummy GSL solver.
    Parameters parameters;

    gsl_multiroot_function testFunction = { &computeAtomTestFunction, 3, &parameters };

    const gsl_multiroot_fsolver_type* solverType = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc( solverType, 3 );
    gsl_multiroot_fsolver_set( solver, &testFunction, initialGuess );

    // Set expected output string.
    std::ostringstream tableRow;
    tableRow << "0  0.1            0.2            0.3            "
             << "1.2            2.3            3.4            "
             << std::endl;

    REQUIRE( printAtomSolverState( 0, solver ) == tableRow.str( ) );
}

TEST_CASE( "Print element", "[print]" )
{
    REQUIRE( printElement( "test", 10 ) == "test      " );    
    REQUIRE( printElement( 1.2345, 10, '.' ) == "1.2345...." );
}

} // namespace tests
} // namespace atom
