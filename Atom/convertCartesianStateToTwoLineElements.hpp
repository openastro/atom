/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H
#define ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H

#include <libsgp4/DateTime.h>
#include <libsgp4/Globals.h>
#include <libsgp4/Tle.h>

namespace atom
{

//! Convert Cartesian state to TLE (Two Line Elements).
/*!
 * Converts a given Cartesian state (position, velocity) to an equivalent TLE. 
 *
 * This function makes use of a root-finder to solve a non-linear system. Locating the root of the 
 * non-linear system corresponds finding the TLE that, when evaluated at its epoch using the 
 * SGP4/SDP4 propagator (Vallado, 2012), yields the target Cartesian state (within tolerance).
 *
 * Details of the underlying non-linear system and algorithm are catalogued by 
 * Kumar, et al. (2014).
 *
 * @sa      evaluateCartesianToTwoLineElementsSystem, DateTime
 * @param   cartesianState Cartesian state [km; km/s]
 * @param   epoch Epoch associated with Cartesian state, stored in a DateTime object
 * @param   solverStatusSummary Status of non-linear solver printed as a table [default: empty] 
 * @param   referenceTle Reference Two Line Elements. This is used as reference to construct the 
 *          effective TLE for the given Cartesian state [default: 0-TLE].
 * @param   earthGravitationalParameter Earth gravitational parameter [km^3 s^-2] [default: mu_SGP]
 * @param   absoluteTolerance Absolute tolerance used to check if root-finder has converged 
 *          [default: 1.0e-5] (see Kumar, et al. (2014) for details on how convergence is tested)
 * @param   relativeTolerance Relative tolerance used to check if root-finder has converged
 *          [default: 1.0e-5] (see Kumar, et al. (2014) for details on how convergence is tested)
 * @param   maximumIterations Maximum number of solver iterations permitted. Once the solver 
 *          reaches this limit, the loop will be broken and the solver status will report that it 
 *          has not converged [default: 100].
 * @return  TLE object that generates target Cartesian state when propagated with SGP4 propagator 
 *          to target epoch.
 */
template< typename Real, typename State >
const Tle convertCartesianStateToTwoLineElements( 
  const State cartesianState,
  const DateTime epoch,
  std::string& solverStatusSummary = "",
  const Tle referenceTle = Tle( ),
  const Real earthGravitationalParameter = kMU,
  const Real absoluteTolerance = 1.0e-5,
  const Real relativeTolerance = 1.0e-5,
  const int maximumIterations = 100 )
{

}

} // namespace atom

#endif // ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H
