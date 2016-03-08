/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Abhishek Agrawal, Delft University of Technology (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <string>
#include <utility>
#include <vector>
#include <stdexcept>

#include <catch.hpp>

#include <libsgp4/DateTime.h>

#include "Atom/atom.hpp"

namespace atom
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector3;
typedef std::pair< Vector3, Vector3 > Velocities;

TEST_CASE( "Execute Atom Solver", "[atom-solver]" )
{
    // Set departure position [km].
    Vector3 departurePosition( 3 );
    departurePosition[ 0 ] = 385.246482619633;
    departurePosition[ 1 ] = 2243.00298468343;
    departurePosition[ 2 ] = 6723.68851012157;

    // Set departure velocity [km/s].
    Vector3 departureVelocity( 3 );
    departureVelocity[ 0 ] = 2.38102833809095;
    departureVelocity[ 1 ] = -6.83464321936554;
    departureVelocity[ 2 ] = 2.02594519285685;

    // Set arrival position [km].
    Vector3 arrivalPosition( 3 );
    arrivalPosition[ 0 ] = -92.709182054428;
    arrivalPosition[ 1 ] = -3085.26470463394;
    arrivalPosition[ 2 ] = -6491.00430424393;

    // Set arrival velocity [km/s].
    Vector3 arrivalVelocity( 3 );
    arrivalVelocity[ 0 ] = -2.39479831205517;
    arrivalVelocity[ 1 ] = 6.3102097995122;
    arrivalVelocity[ 2 ] = -3.09055920479061;

    // Set departure epoch.
    DateTime departureEpoch( 2016, 1, 11, 4, 53, 22 );

    // Time-of-flight [s].
    const Real timeOfFlight = 2830.0;

    SECTION( "Test case with no iterations" )
    {
        // Execute Atom solver to compute departure and arrival velocities that bridge the given
        // positions.
        std::string dummyString = "";
        int numberOfIterations = 0;
        const Velocities velocities = executeAtomSolver( departurePosition,
                                                         departureEpoch,
                                                         arrivalPosition,
                                                         timeOfFlight,
                                                         departureVelocity,
                                                         dummyString,
                                                         numberOfIterations );

        // Check that departure and arrival velocities match results.
        for ( int i = 0; i < 3; i++ )
        {
            REQUIRE( departureVelocity[ i ] == Approx( velocities.first[ i ] ).epsilon( 1.0e-6 ) );
            REQUIRE( arrivalVelocity[ i ] == Approx( velocities.second[ i ] ).epsilon( 1.0e-6 ) );
        }

        // Check that no iterations are required.
        REQUIRE( numberOfIterations == 0 );
    }

    SECTION( "Test arbitrary case" )
    {
        // Set initial guess for departure velocity [km/s]. 
        Vector3 departureVelocityGuess( 3 );
        departureVelocityGuess[ 0 ] = 2.39618773382003;
        departureVelocityGuess[ 1 ] = -6.83412563589514;
        departureVelocityGuess[ 2 ] = 2.03013574194101;

        // Execute Atom solver to compute departure and arrival velocities that bridge the given
        // positions.
        std::string dummyString = "";
        int numberOfIterations = 0;
        const Velocities velocities = executeAtomSolver( departurePosition,
                                                         departureEpoch,
                                                         arrivalPosition,
                                                         timeOfFlight,
                                                         departureVelocityGuess,
                                                         dummyString,
                                                         numberOfIterations );

        // Check that departure and arrival velocities match results.
        for ( int i = 0; i < 3; i++ )
        {
            REQUIRE( departureVelocity[ i ] == Approx( velocities.first[ i ] ).epsilon( 1.0e-6 ) );
            REQUIRE( arrivalVelocity[ i ] == Approx( velocities.second[ i ] ).epsilon( 1.0e-6 ) );
        }

        // Check that no iterations are required.
        REQUIRE( numberOfIterations == 2 );
    }
}

// This is a legacy test which should result in an exception being thrown since the tsince bug was fixed
TEST_CASE( "Execute Legacy Atom solver Test", "[atom-solver]")
{
    // Set departure position [km].
    Vector3 departurePosition( 3 );
    departurePosition[ 0 ] = -3680.20448307549;
    departurePosition[ 1 ] = -2573.44661796266;
    departurePosition[ 2 ] = 5800.72628190982;

    // Set departure velocity [km/s].
    Vector3 departureVelocity( 3 );
    departureVelocity[ 0 ] = 6.44661660560979;
    departureVelocity[ 1 ] = -1.14788435945363;
    departureVelocity[ 2 ] = 3.44659369332744;

    // Set arrival position [km].
    Vector3 arrivalPosition( 3 );
    arrivalPosition[ 0 ] = 4496.59209320659;
    arrivalPosition[ 1 ] = 2339.99159226651;
    arrivalPosition[ 2 ] = -5455.56445926525;

    // Set arrival velocity [km/s].
    Vector3 arrivalVelocity( 3 );
    arrivalVelocity[ 0 ] = -5.7447123573464;
    arrivalVelocity[ 1 ] = 1.63941146365299;
    arrivalVelocity[ 2 ] = -4.17792707643158;

    // Set departure epoch.
    DateTime departureEpoch( 63548650522376360 );

    // Time-of-flight [s].
    const Real timeOfFlight = 1000.0;

    SECTION( "Test case with no iterations" )
    {
        // Execute Atom solver to compute departure and arrival velocities that bridge the given
        // positions.
        std::string dummyString = "";
        int numberOfIterations = 0;
        
        // check that an exception is thrown
        REQUIRE_THROWS( executeAtomSolver( departurePosition,
                                           departureEpoch,
                                           arrivalPosition,
                                           timeOfFlight,
                                           departureVelocity,
                                           dummyString,
                                           numberOfIterations ) );
    }

    SECTION( "Test arbitrary case" )
    {
        // Set initial guess for departure velocity [km/s]. Arbitrary values are added to the
        // expected departure velocity.
        Vector3 departureVelocityGuess( 3 );
        departureVelocityGuess[ 0 ] = departureVelocity[ 0 ] + 0.013;
        departureVelocityGuess[ 1 ] = departureVelocity[ 1 ] - 0.074;
        departureVelocityGuess[ 2 ] = departureVelocity[ 2 ] + 0.026;
        
        // Execute Atom solver to compute departure and arrival velocities that bridge the given
        // positions.
        std::string dummyString = "";
        int numberOfIterations = 0;

        // check that an exception is thrown
        REQUIRE_THROWS( executeAtomSolver( departurePosition,
                                           departureEpoch,
                                           arrivalPosition,
                                           timeOfFlight,
                                           departureVelocity,
                                           dummyString,
                                           numberOfIterations ) );
    }
}

} // namespace tests
} // namespace atom
