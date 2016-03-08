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

// This is a legacy test which should result in an exception being thrown since the tsince bug was fixed
TEST_CASE( "Execute Atom solver", "[atom-solver]")
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
