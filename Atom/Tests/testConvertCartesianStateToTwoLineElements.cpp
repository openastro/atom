/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#define INTEGER int
#define REAL double

#include <string>
#include <vector>

#include <catch.hpp>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>
 
#include "Atom/convertCartesianStateToTwoLineElements.hpp"

namespace atom
{
namespace tests
{

typedef std::vector< REAL > Vector;

TEST_CASE( "Convert Cartesian state to Two-Line-Elements", "[cartesian-to-TLE]")
{
    // Set target cartesian state [km; km/s].
    Vector cartesianState( 6 );
    cartesianState[ 0 ] = -7.1e3;
    cartesianState[ 1 ] = 2.7e3;
    cartesianState[ 2 ] = 1.3e3;
    cartesianState[ 3 ] = -2.5;
    cartesianState[ 4 ] = -5.5;
    cartesianState[ 5 ] = 5.5;

    // Convert Cartesian state to TLE. Note that the epoch is arbitrary.
    Tle convertedTle = convertCartesianStateToTwoLineElements< INTEGER, REAL, Vector >( 
        cartesianState, DateTime( ) );

    // Propagate the converted TLE to the epoch of the TLE. This generates a Cartesian state.
    SGP4 sgp4( convertedTle );
    Eci recomputedCartesianState = sgp4.FindPosition( 0.0 );

    // Check if recomputed Cartesian state matches the target state.
    REQUIRE( recomputedCartesianState.Position( ).x 
             == Approx( cartesianState[ 0 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedCartesianState.Position( ).y 
             == Approx( cartesianState[ 1 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedCartesianState.Position( ).z 
             == Approx( cartesianState[ 2 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedCartesianState.Velocity( ).x 
             == Approx( cartesianState[ 3 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedCartesianState.Velocity( ).y 
             == Approx( cartesianState[ 4 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedCartesianState.Velocity( ).z 
             == Approx( cartesianState[ 5 ] ).epsilon( 1.0e-8 ) );
}

} // namespace tests
} // namespace atom
