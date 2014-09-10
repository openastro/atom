/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <utility>

#include <Eigen/Core>

#include <libsgp4/DateTime.h>
#include <libsgp4/Globals.h>
#include <libsgp4/Tle.h>

#include "Atom/generalDefinitions.h"
 
namespace atom
{

//! Execute Atom solver.
/*!
 * Executes Atom solver to find the transfer orbit connecting two positions, at two specified 
 * epochs. 
 *
 * The Atom solver is an analog of the Lambert solver (Lancaster and Blanchard, 1969;
 * Gooding, 1990; Izzo, 2014), that aims to find the conic section that bridges two positions, at
 * given epochs, by using impulsive manveuvers (DeltaV maneuvers) at departure and arrival. The
 * Atom solver aims to solver a similar orbital transfer, subject to perturbations. The 
 * perturbations taken into account are those encoded in the SGP4/SDP4 propagators (Vallado, 2006).
 *
 * Since the Atom solver makes use fo the SGP4/SDP4 propagators, it can currently only solve for
 * perturbed transfers around the Earth. As a result, the Earth's gravitational parameter is fixed,
 * as specified by the SGP4/SDP4 propagators (Vallado, 2006).
 * 
 * \param departureState Cartesian state vector at departure point. Note that this is BEFORE any
 *          departure Delta V is applied! [m; m/s].
 * \param departureEpoch Modified Julian Date (MJD) of departure.
 * \param arrivalState Cartesian state vector at arrival point. Note that this is the target 
 *          arrival state AFTER any arrival Delta V is applied! [m; m/s].
 * \param timeOfFlight Time-of-flight (TOF) for orbital transfer [s].
 * \param departureDeltaVInitialGuess Initial guess for the departure Delta V. This serves as
 *          initial guess for the internal root-finding procedure [m/s].
 * \param solverStatusSummary Summary table containing status of non-linear solver per iteration. 
 * \param referenceDepartureTle Reference TLE at departure. This is used as reference when 
 *          converting the Cartesian state at departure to a TLE. [default: 0-TLE].
 * \param earthGravitationalParameter Earth gravitational parameter [m^3 s^-2] [default: mu_SGP].  
 * \param convergenceTolerance Convergence tolerance for non-linear solver [default: 1.0e-8].
 * \param maximumIterations Maximum number of iterations permitted for non-linear solver 
 *          [default: 100].
 * \return Departure and arrival DeltaVs (stored in that order).
 */
const std::pair< Eigen::Vector3d, Eigen::Vector3d > executeAtomSolver( 
    const Vector6d departureState, 
    const DateTime departureEpoch,
    const Vector6d arrivalState, 
    const double timeOfFlight, 
    const Eigen::Vector3d departureDeltaVInitialGuess,
    std::string& solverStatusSummary = emptyString,        
    const Tle referenceDepartureTle = Tle( ),
    const double earthGravitationalParameter = kMU * 1.0e9,    
    const double convergenceTolerance = 1.0e-8,
    const int maximumIterations = 100 );

} // namespace atom

#endif // ATOM_H
