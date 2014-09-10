/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>  
#include <libsgp4/SGP4.h>  
#include <libsgp4/Tle.h>     

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>  

#include "Atom/atom.h"

//! Execute main-function.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Declare using-statements.
    using std::string;

    using namespace tudat::basic_astrodynamics::unit_conversions;

    using namespace atom;

    //////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Departure epoch.
    DateTime departureEpoch( 2014, 5, 23, 15, 0, 0 );

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

    //////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////

    // Compute derived parameters.

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

    // Convert Cartesian departure state to a new TLE object.
    string solverSummaryTable = "";

    // const Tle newDepartureTle = convertCartesianStateToTwoLineElements( 
    //     departureState, 
    //     tleDepartureObject.Epoch( ),
    //     solverSummaryTable,
    //     Tle( ), 
    //     earthGravitationalParameter );

// executeAtomSolver( 
//      departureState, 
//      tleDepartureObject.Epoch( ),
//      departureState, 
//      tleDepartureObject.Epoch( ),
//      1.0, 
//      Eigen::Vector3d( 0.0, 10.0, 0.0 ),
//      solverSummaryTable );    

std::cout << solverSummaryTable << std::endl; 

// std::cout << departureState << std::endl;
// std::cout << std::endl;
// std::cout << tleDepartureObject << std::endl;     
// std::cout << std::endl;   
// std::cout << newDepartureTle << std::endl;        

    //////////////////////////////////////////////////////////////////////////////////////////////

	return EXIT_SUCCESS;
}