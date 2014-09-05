/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include "Atom/atom.h"

namespace atom
{

//! Execute Atom solver.
const std::pair< Eigen::Vector3d, Eigen::Vector3d > executeAtomSolver( 
    const Vector6d departureState, 
    const double departureEpoch,
    const Vector6d arrivalState, 
    const double arrivalEpoch,
    const double timeOfFlight, 
    const Tle referenceDepartureTle,
    const Tle referenceArrivalTle,
    std::string& solverStatusSummary,    
    const double convergenceTolerance,
    const int maximumIterations )
{

    // TEMP
    return std::make_pair< Eigen::Vector3d, Eigen::Vector3d >( Eigen::Vector3d( 0.0, 0.0, 0.0 ), Eigen::Vector3d( 0.0, 0.0, 0.0 ) );
}

} // namespace atom
