/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE Atom

#include <string>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>  

#include "Atom/atom.h"
#include "Atom/generalDefinitions.h"

namespace atom
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_atom )

//! Test implementation of Atom solver.
BOOST_AUTO_TEST_CASE( testAtomSolver )
{
    std::cout << "Atom test" << std::endl;

    // Declare using-statements.
    using std::string;
    using namespace tudat::basic_astrodynamics::unit_conversions;

    // Set up TLE strings.
    const string nameDepartureObject = "TEST";
    const string line1DepartureObject
        = "1 00341U 62029B   12295.78875316 -.00000035  00000-0  68936-4 0  4542";
    const string line2DepartureObject
        = "2 00341 044.7940 159.4482 2422241 038.6572 336.4786 09.14200419679435"; 

    // Time-of-flight [s].
    const double timeOfFlight = 1000.0;

    // Departure Delta V [m/s].
    const Eigen::Vector3d departureDeltaV( 10.0, 0.0, 0.0 );    

    // Set up TLE objects.
    const Tle tleDepartureObject( 
        nameDepartureObject, line1DepartureObject, line2DepartureObject );

    // Set up SGP4 propagator for TLE objects.
    const SGP4 sgp4DepartureObject( tleDepartureObject );

    // Propagate object to epoch of TLE.
    Eci sgp4DepartureState = sgp4DepartureObject.FindPosition( 0.0 );    

    // Store propagated state at TLE epoch as Cartesian state.
    const Vector6d departureState 
        = ( Eigen::VectorXd( 6 ) 
                << convertKilometersToMeters( sgp4DepartureState.Position( ).x ),
                   convertKilometersToMeters( sgp4DepartureState.Position( ).y ),
                   convertKilometersToMeters( sgp4DepartureState.Position( ).z ),
                   convertKilometersToMeters( sgp4DepartureState.Velocity( ).x ),
                   convertKilometersToMeters( sgp4DepartureState.Velocity( ).y ),
                   convertKilometersToMeters( sgp4DepartureState.Velocity( ).z ) ).finished( );

    // Propagate object to arrival epoch.
    Eci sgp4ArrivalState = sgp4DepartureObject.FindPosition( timeOfFlight );    

    // Store propagated state at arrival epoch as Cartesian state.
    const Vector6d arrivalState 
        = ( Eigen::VectorXd( 6 ) 
                << convertKilometersToMeters( sgp4ArrivalState.Position( ).x ),
                   convertKilometersToMeters( sgp4ArrivalState.Position( ).y ),
                   convertKilometersToMeters( sgp4ArrivalState.Position( ).z ),
                   convertKilometersToMeters( sgp4ArrivalState.Velocity( ).x ),
                   convertKilometersToMeters( sgp4ArrivalState.Velocity( ).y ),
                   convertKilometersToMeters( sgp4ArrivalState.Velocity( ).z ) ).finished( );

    const DeltaVs deltaVs = executeAtomSolver( departureState, 
                                               tleDepartureObject.Epoch( ),
                                               arrivalState, 
                                               timeOfFlight, 
                                               Eigen::Vector3d( 0.0, 0.0, 0.0 ) ); 

    std::cout << "Departure Delta V: " << deltaVs.first << std::endl;   
    std::cout << "Arrival Delta V: " << deltaVs.second << std::endl;   

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace atom
