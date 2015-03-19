/*
 * Copyright (c) 2014-2015 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <string>

#include <catch.hpp>

#include <Eigen/Core>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include "Atom/convertCartesianStateToTwoLineElements.hpp"

namespace atom
{
namespace tests
{

typedef double Real;
typedef Eigen::Matrix< Real, Eigen::Dynamic, 1 > Vector;

TEST_CASE( "Convert Cartesian state to Two-Line-Elements", "[cartesian-to-TLE]")
{
    // Set target cartesian state [km; km/s].
    Vector cartesianElements( 6 );
    cartesianElements[ 0 ] = -7.1e3;
    cartesianElements[ 1 ] = 2.7e3;
    cartesianElements[ 2 ] = 1.3e3;
    cartesianElements[ 3 ] = -2.5;
    cartesianElements[ 4 ] = -5.5;
    cartesianElements[ 5 ] = 5.5;

    // Convert Cartesian state to TLE. Note that the epoch is arbitrary.
    const Tle convertedTle = convertCartesianStateToTwoLineElements< Real, Vector >(
        cartesianElements, DateTime( ) );

    // Propagate the converted TLE to the epoch of the TLE. This generates a Cartesian state.
    SGP4 sgp4( convertedTle );
    Eci recomputedcartesianElements = sgp4.FindPosition( 0.0 );

    // Check if recomputed Cartesian state matches the target state.
    REQUIRE( recomputedcartesianElements.Position( ).x
             == Approx( cartesianElements[ 0 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedcartesianElements.Position( ).y
             == Approx( cartesianElements[ 1 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedcartesianElements.Position( ).z
             == Approx( cartesianElements[ 2 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedcartesianElements.Velocity( ).x
             == Approx( cartesianElements[ 3 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedcartesianElements.Velocity( ).y
             == Approx( cartesianElements[ 4 ] ).epsilon( 1.0e-8 ) );
    REQUIRE( recomputedcartesianElements.Velocity( ).z
             == Approx( cartesianElements[ 5 ] ).epsilon( 1.0e-8 ) );
}

} // namespace tests
} // namespace atom
