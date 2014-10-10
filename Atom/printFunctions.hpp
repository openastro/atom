/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <gsl/gsl_multiroots.h>

#ifndef ATOM_PRINT_FUNCTIONS_H
#define ATOM_PRINT_FUNCTIONS_H

namespace atom
{

//! Print solver summary table header.
/*!
 * Prints header to string for table containing summary of status of non-linear solver used to 
 * convert a Cartesian state to a TLE.
 *
 * @sa convertCartesianStateToTwoLineElements 
 * @return String containing table header for non-linear solver status.
 */
inline std::string printSolverStateTableHeader( );

//! Print summary of current state of non-linear solver.
/*!
 * Prints current state of non-linear solver used to convert a Cartesian state to a TLE, as row 
 * for a summary table.
 *
 * @sa convertCartesianStateToTwoLineElements
 * @param  iteration Current iteration of solver
 * @param  solver    Pointer to GSL solver
 * @return           String containing row-data for non-linear solver status summary table
 */
inline std::string printSolverState( const int iteration, gsl_multiroot_fsolver* solver );

//! Print data element to console.
/*!
 * Prints a specified data element to string, given a specified width and a separator character.
 * This function is auxilliary to the print-functions used to the print the state of the non-linear
 * solver.
 * 
 * @sa printSolverStateTableHeader, printSolverState
 * @tparam DataType  Type for specified data element
 * @param  datum     Specified data element to print
 * @param  width     Width of datum printed to console, in terms of number of characters
 * @param  separator Separator character, e.g., ","
 * @return           String containing printed data element
 */
template< typename DataType > 
inline std::string printElement( const DataType datum, const int width, const char separator );

//! Print solver summary table header.
inline std::string printSolverStateTableHeader( )
{
    std::ostringstream headerBuffer;
    headerBuffer << printElement( "#", 3, ' ' )
                 << printElement( "a", 15, ' ' )
                 << printElement( "e", 15, ' ' )
                 << printElement( "i", 15, ' ' )
                 << printElement( "AoP", 15, ' ' )
                 << printElement( "RAAN", 15, ' ' )
                 << printElement( "TA", 15, ' ' )
                 << printElement( "f1", 15, ' ' )
                 << printElement( "f2", 15, ' ' )
                 << printElement( "f3", 15, ' ' )
                 << printElement( "f4", 15, ' ' )
                 << printElement( "f5", 15, ' ' )
                 << printElement( "f6", 15, ' ' )
                 << std::endl;
    return headerBuffer.str( );    
}

//! Print summary of current state of non-linear solver.
/*!
 * Prints current state of non-linear solver used to convert a Cartesian state to a TLE, as row 
 * for a summary table.
 *
 * @sa convertCartesianStateToTwoLineElements
 * @param  iteration Current iteration of solver
 * @param  solver    Pointer to GSL solver
 * @return           String containing row-data for non-linear solver status summary table
 */
inline std::string printSolverState( const int iteration, gsl_multiroot_fsolver* solver )
{
    std::ostringstream buffer;
    buffer << printElement( iteration, 3, ' ' )
           << printElement( gsl_vector_get( solver->x, 0 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 1 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 2 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 3 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 4 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 5 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 0 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 1 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 2 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 3 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 4 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 5 ), 15, ' ' )
           << std::endl;    
    return buffer.str( );
}

//! Print data element to console.
template< typename DataType > 
inline std::string printElement( const DataType datum, const int width, const char separator )
{
    std::ostringstream buffer;
    buffer << std::left << std::setw( width ) << std::setfill( separator ) << datum;
    return buffer.str( );
}

} // namespace atom

#endif // ATOM_PRINT_FUNCTIONS_H
