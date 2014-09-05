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
 * \param departureState Cartesian state vector at departure point [m; m/s].
 * \param departureEpoch Modified Julian Date (MJD) of departure.
 * \param arrivalState Cartesian state vector at arrival point [m; m/s].
 * \param arrivalEpoch Modified Julian Date (MJD) of arrival.
 * \param timeOfFlight Time-of-flight (TOF) for orbital transfer [s].
 * \param referenceDepartureTle Reference TLE at departure. This is used as reference when 
 *          converting the Cartesian state at departure to a TLE. [default: 0-TLE].
 * \param referenceArrivalTle Reference TLE at arrival. This is used as reference when 
 *          converting the Cartesian state at arrival to a TLE. [default: 0-TLE].
 * \param solverStatusSummary Summary table containing status of non-linear solver per iteration. 
 * \param convergenceTolerance Convergence tolerance for non-linear solver [default: 1.0e-8].
 * \param maximumIterations Maximum number of iterations permitted for non-linear solver 
 *          [default: 100].
 * \return Departure and arrival DeltaVs (stored in that order).s
 */
const std::pair< Eigen::Vector3d, Eigen::Vector3d > executeAtomSolver( 
    const Vector6d departureState, 
    const double departureEpoch,
    const Vector6d arrivalState, 
    const double arrivalEpoch,
    const double timeOfFlight, 
    const Tle referenceDepartureTle = Tle( ),
    const Tle referenceArrivalTle = Tle( ),
    std::string& solverStatusSummary = emptyString,    
    const double convergenceTolerance = 1.0e-8,
    const int maximumIterations = 100 );

} // namespace atom

#endif // ATOM_H
