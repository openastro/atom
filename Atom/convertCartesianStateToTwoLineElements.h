/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H
#define ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
 
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include <libsgp4/Globals.h>
#include <libsgp4/Tle.h>

#include "Atom/generalDefinitions.H"

namespace atom
{

//! Convert Cartesian state to TLE (Two Line Elements).
/*!
 * Converts a given Cartesian state (position, velocity) to an equivalent TLE. 
 *
 * This function makes use of a root-finder to solve a non-linear system. Locating the root of the 
 * non-linear system corresponds finding the TLE that, when evaluated at its epoch using the 
 * SGP4/SDP4 propagator, yields the target Cartesian state (within tolerance).
 *
 * Details of the underlying non-linear system and algorithm are catalogued by 
 * Kumar, et al. (2014).
 *
 * \param cartesianState Cartesian state.
 * \param epoch Epoch associated with Cartesian state.
 * \param earthGravitationalParameter Earth gravitational parameter [m^3 s^-2] [default: mu_SGP].
 * \param referenceTle Reference Two Line Elements. This is used as reference to construct the 
 *         effective TLE for the given Cartesian state [default: 0-TLE].
 * \param solverStatusSummary Status of the non-linear solver printed as a table [default: empty].
 * \param tolerance Tolerance used to check if root-finder has converged [default: 1.0e-8].
 * \param maximumIterations Maximum number of solver iterations permitted. Once the solver reaches
 *          this limit, the loop will be broken and the solver status will report that it has not
 *          converged [default: 100]. 
 * \return TLE object that generates target Cartesian state when propagated with SGP4 propagator to
 *           target epoch.
 * \sa evaluateCartesianToTwoLineElementsSystem()
 */
const Tle convertCartesianStateToTwoLineElements( 
  const Vector6d cartesianState,
  const DateTime epoch,
  std::string& solverStatusSummary = emptyString,
  const Tle referenceTle = Tle( ),
  const double earthGravitationalParameter = kMU * 1.0e9,
  const double tolerance = 1.0e-6,
  const int maximumIterations = 100 );

//! Evaluate system of non-linear equations for converting Cartesian state to TLE.
/*!
 * Evaluates system of non-linear equations to find TLE corresponding with target Cartesian 
 * state. The system of non-linear equations used is:
 *  \f[ 
 *      F = 0 = \bar{x}^{new} - \bar{x}^{target}
 *  \f]
 * where \f$\bar{x}^{new}\f$ is the new Cartesian state computed by updating the TLE mean elements 
 * and propagating the TLE using the SGP4 propagator, and \f$\bar{x}^{target}\f$ is the target
 * Cartesian state.
 * \param independentVariables Vector of independent variables used by the root-finder.
 * \param parameters Parameters required to compute the objective function.
 * \param functionValues Vector of computed function values. 
 */
int evaluateCartesianToTwoLineElementsSystem( const gsl_vector* independentVariables, 
                                              void* parameters, 
                                              gsl_vector* functionValues );

//! Update TLE mean elements.
/*!
 * Updates mean elements stored in TLE based on current osculating elements. This function 
 * uses the osculating elements to replace mean elements (converts units and computes mean anomaly
 * and mean motion).
 * \param newKeplerianElements New Keplerian elements.
 * \param oldTle TLE in which the mean elements are to be replaced.
 * \param earthGravitationalParameter Earth gravitational parameter [m^3 s^-2].
 * \return New TLE with mean elements updated.
 */
const Tle updateTleMeanElements( const Vector6d newKeplerianElements, 
    							               const Tle oldTle,
    							               const double earthGravitationalParameter );

//! Container of parameters used by Cartesian-to-TLE objective function.
struct CartesianToTwoLineElementsObjectiveParameters
{ 
public:

    //! Constructor taking parameter values.
    CartesianToTwoLineElementsObjectiveParameters( 
        const Vector6d aTargetState,
        const double anEarthGravitationalParameter,
        const Tle someReferenceTle )
        : targetState( aTargetState ),
          earthGravitationalParameter( anEarthGravitationalParameter ),
 		  referenceTle( someReferenceTle )
    { }

    //! Target state in Cartesian elements.
    const Eigen::VectorXd targetState;

    //! Earth gravitational parameter [kg m^3 s^-2].
    const double earthGravitationalParameter;

    //! Reference TLE.
    const Tle referenceTle;

protected:

private:
};

//! Print header for table containing summary of non-linear solver state.
/*!
 * Prints header to string for table containing summary of status of non-linear solver used to 
 * convert a Cartesian state to a TLE.
 * \return String containing table header for non-linear solver status.
 * \sa convertCartesianStateToTwoLineElements()
 */
std::string printSolverStateTableHeader( );

//! Print current state of non-linear solver for summary table.
/*!
 * Prints current state of non-linear solver used to convert a Cartesian state to a TLE, as row for 
 * a summary table.
 * \param iteration Current iteration of solver.
 * \param solver Pointer to GSL solver.
 * \return String containing row-data for non-linear solver status summary table.
 */
std::string printSolverState( int iteration, gsl_multiroot_fsolver* solver );

//! Print data element to console.
/*!
 * Prints a specified data element to string, given a specified width and a separator character.
 * This function is auxilliary to the print-functions used to the print the state of the non-linear
 * solver.
 * \tparam DataType Type for specified data element.
 * \param datum Specified data element to print.
 * \param width Width of datum printed to console, in terms of number of characters.
 * \param separator Separator character, e.g., ",".
 * \return String containing printed data element.
 * \sa printSolverState().
 */
template< typename DataType > 
inline std::string printElement( const DataType datum, const int width, const char separator )
{
    std::ostringstream buffer;
    buffer << std::left << std::setw( width ) << std::setfill( separator ) << datum;
    return buffer.str( );
}

} // namespace atom

#endif // ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H
