/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#define INTEGER int
#define REAL double

#include <iomanip>
#include <limits>

#include <string>
#include <utility>
#include <vector>

#include <catch.hpp>

#include <libsgp4/DateTime.h>
 
#include "Atom/atom.hpp"

namespace atom
{
namespace tests
{

typedef std::vector< REAL > Vector;
typedef std::pair< Vector, Vector > Velocities;

TEST_CASE( "Execute Atom solver", "[atom-solver]")
{
    // Set departure position [km].
    Vector departurePosition( 3 );
    departurePosition[ 0 ] = -3680.20448307549;
    departurePosition[ 1 ] = -2573.44661796266;
    departurePosition[ 2 ] = 5800.72628190982;

    // Set departure velocity [km/s].
    Vector departureVelocity( 3 );
    departureVelocity[ 0 ] = 6.44661660560979;
    departureVelocity[ 1 ] = -1.14788435945363;
    departureVelocity[ 2 ] = 3.44659369332744;

    // Set arrival position [km].
    Vector arrivalPosition( 3 );
    arrivalPosition[ 0 ] = 4496.59209320659;
    arrivalPosition[ 1 ] = 2339.99159226651;
    arrivalPosition[ 2 ] = -5455.56445926525;

    // Set arrival velocity [km/s].
    Vector arrivalVelocity( 3 );
    arrivalVelocity[ 0 ] = -5.7447123573464;
    arrivalVelocity[ 1 ] = 1.63941146365299;
    arrivalVelocity[ 2 ] = -4.17792707643158;

    // Set departure epoch.
    DateTime departureEpoch( 63548650522376360 );

    // Time-of-flight [s].
    const REAL timeOfFlight = 1000.0;
            
    SECTION( "Test case with no iterations" )
    {
        // Execute Atom solver to compute departure and arrival velocities that bridge the given
        // positions.
        std::string dummyString = "";
        INTEGER numberOfIterations = 0;
        const Velocities velocities = executeAtomSolver< INTEGER, REAL, Vector >( 
            departurePosition, 
            departureEpoch, 
            arrivalPosition, 
            timeOfFlight, 
            departureVelocity,
            dummyString,
            numberOfIterations ); 

        // Check that departure and arrival velocities match results.
        for ( INTEGER i = 0; i < 3; i++ )
        {
            REQUIRE( departureVelocity[ i ] == Approx( velocities.first[ i ] ).epsilon( 1.0e-6 ) );
            REQUIRE( arrivalVelocity[ i ] == Approx( velocities.second[ i ] ).epsilon( 1.0e-6 ) );
        }

        // Check that no iterations are required.
        REQUIRE( numberOfIterations == 0 );
    }                                       
}

} // namespace tests
} // namespace atom
