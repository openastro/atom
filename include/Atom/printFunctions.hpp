/*    
 * Copyright (c) 2014 K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
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

//! Print Cartesian-state-to-TLE converter solver summary table header.
/*!
 * Prints header to string for table containing summary of status of non-linear solver used to 
 * convert a Cartesian state to a TLE.
 *
 * @sa convertCartesianStateToTwoLineElements 
 * @return String containing table header for non-linear solver status.
 */
inline std::string printCartesianToTleSolverStateTableHeader( );

//! Print summary of current state of non-linear solver for Cartesian-state-to-TLE converter.
/*!
 * Prints current state of non-linear solver used to convert a Cartesian state to a TLE, as row 
 * for a summary table.
 *
 * @sa convertCartesianStateToTwoLineElements
 * @param  iteration Current iteration of solver
 * @param  solver    Pointer to GSL solver
 * @return           String containing row-data for non-linear solver status summary table
 */
inline std::string printCartesianToTleSolverState( 
  const int iteration, gsl_multiroot_fsolver* solver );

//! Print Atom solver summary table header.
/*!
 * Prints header to string for table containing summary of status of non-linear solver used to 
 * execute Atom solver.
 *
 * @sa executeAtomSolver 
 * @return String containing table header for non-linear solver status.
 */
inline std::string printAtomSolverStateTableHeader( );

//! Print summary of current state of non-linear solver for Atom solver.
/*!
 * Prints current state of non-linear solver used to execute Atom solver, as row for a summary
 * table.
 *
 * @sa executeAtomSolver
 * @param  iteration Current iteration of solver
 * @param  solver    Pointer to GSL solver
 * @return           String containing row-data for non-linear solver status summary table
 */
inline std::string printAtomSolverState( 
  const int iteration, gsl_multiroot_fsolver* solver );

//! Print data element to console.
/*!
 * Prints a specified data element to string, given a specified width and a filler character.
 * This function is auxilliary to the print-functions used to the print the state of the non-linear
 * solver.
 * 
 * @sa printSolverStateTableHeader, printSolverState
 * @tparam DataType  Type for specified data element
 * @param  datum     Specified data element to print
 * @param  width     Width of datum printed to console, in terms of number of characters
 * @param  filler    Character used to fill fixed-width, [default: ' ']
 * @return           String containing printed data element
 */
template< typename DataType > 
inline std::string printElement( const DataType datum, const int width, const char filler = ' ' );

//! Print Cartesian-state-to-TLE converter solver summary table header.
inline std::string printCartesianToTleSolverStateTableHeader( )
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

//! Print summary of current state of non-linear solver for Cartesian-state-to-TLE converter.
inline std::string printCartesianToTleSolverState( 
  const int iteration, gsl_multiroot_fsolver* solver )
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

//! Print Atom solver summary table header.
inline std::string printAtomSolverStateTableHeader( )
{
    std::ostringstream headerBuffer;
    headerBuffer << printElement( "#", 3, ' ' )
                 << printElement( "v1_x", 15, ' ' )
                 << printElement( "v1_y", 15, ' ' )
                 << printElement( "v1_z", 15, ' ' )
                 << printElement( "f1", 15, ' ' )
                 << printElement( "f2", 15, ' ' )
                 << printElement( "f3", 15, ' ' )
                 << std::endl;
    return headerBuffer.str( );   
}

//! Print summary of current state of non-linear solver for Atom solver.
inline std::string printAtomSolverState( 
  const int iteration, gsl_multiroot_fsolver* solver )
{
    std::ostringstream buffer;
    buffer << printElement( iteration, 3, ' ' )
           << printElement( gsl_vector_get( solver->x, 0 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 1 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 2 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 0 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 1 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 2 ), 15, ' ' )
           << std::endl;    
    return buffer.str( );
}

//! Print data element to console.
template< typename DataType > 
inline std::string printElement( const DataType datum, const int width, const char filler )
{
    std::ostringstream buffer;
    buffer << std::left << std::setw( width ) << std::setfill( filler ) << datum;
    return buffer.str( );
}

} // namespace atom

#endif // ATOM_PRINT_FUNCTIONS_H
