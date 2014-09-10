/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>  

#include "Atom/atom.h"
#include "Atom/convertCartesianStateToTwoLineElements.h"

namespace atom
{

//! Execute Atom solver.
const std::pair< Eigen::Vector3d, Eigen::Vector3d > executeAtomSolver( 
    const Vector6d departureState, 
    const DateTime departureEpoch,
    const Vector6d arrivalState, 
    const double timeOfFlight, 
    const Eigen::Vector3d departureDeltaVInitialGuess,
    std::string& solverStatusSummary,        
    const Tle referenceDepartureTle,
    const double earthGravitationalParameter,    
    const double convergenceTolerance,
    const int maximumIterations )
{
    // Set departure Delta V.
    Eigen::Vector3d departureDeltaV = departureDeltaVInitialGuess;

    // Compute post-maneuver departure state (after Delta V is applied).
    Vector6d postManeuverDepartureState 
        = departureState + ( Eigen::VectorXd( 6 ) 
                             << Eigen::Vector3d( ) , departureDeltaV  ).finished( );

    // Convert post-maneuver departure state to effective TLE.
    const Tle departureTle = convertCartesianStateToTwoLineElements( 
        postManeuverDepartureState, departureEpoch, solverStatusSummary );

    // Propagate departure TLE by time-of-flight using SGP4 propagator.
    SGP4 sgp4( departureTle );
    Eci preManeuverSgp4ArrivalState = sgp4.FindPosition( timeOfFlight );

    // Store propagated state at TLE epoch as Cartesian state.
    using namespace tudat::basic_astrodynamics::unit_conversions;    
    const Vector6d preManeuverArrivalState 
        = ( Eigen::VectorXd( 6 ) 
            << convertKilometersToMeters( preManeuverSgp4ArrivalState.Position( ).x ),
               convertKilometersToMeters( preManeuverSgp4ArrivalState.Position( ).y ),
               convertKilometersToMeters( preManeuverSgp4ArrivalState.Position( ).z ),
               convertKilometersToMeters( preManeuverSgp4ArrivalState.Velocity( ).x ),
               convertKilometersToMeters( preManeuverSgp4ArrivalState.Velocity( ).y ),
               convertKilometersToMeters( 
                preManeuverSgp4ArrivalState.Velocity( ).z ) ).finished( );  

    // Compute Delta V applied at arrival.
    Eigen::Vector3d arrivalDeltaV 
        = preManeuverArrivalState.segment( 3, 3 ) - arrivalState.segment( 3, 3 );        

    // Return departure and arrival Delta Vs.
    return std::make_pair< Eigen::Vector3d, Eigen::Vector3d >( departureDeltaV, arrivalDeltaV );
}

} // namespace atom
