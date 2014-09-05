/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h> 
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>  
#include <TudatCore/Basics/testMacros.h>

#include "Atom/convertCartesianStateToTwoLineElements.h"
#include "Atom/generalDefinitions.h"

namespace atom
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_cartesian_to_tle_conversion )

//! Test implementation of Cartesian-state-to-TLE converter.
BOOST_AUTO_TEST_CASE( testCartesianStateToTleConverter )
{
    // Set Cartesian state.
    Vector6d cartesianState = ( Eigen::VectorXd( 6 ) << -7.1e6, 2.7e6, 1.3e6, 
                                                        -2.5e3, -5.5e3, 5.5e3 ).finished( );

    // Convert Cartesian state to TLE.
    Tle cartesianStateTle 
        = atom::convertCartesianStateToTwoLineElements( cartesianState, DateTime( ) );

    // Use SGP4 propagator to convert TLE back to Cartesian state.
    SGP4 sgp4( cartesianStateTle );
    Eci propagatedState = sgp4.FindPosition( 0.0 );    

    // Store recomputed Cartesian state.
    using namespace tudat::basic_astrodynamics::unit_conversions;    
    const Vector6d recomputedCartesianState 
        = ( Eigen::VectorXd( 6 ) 
            << convertKilometersToMeters( propagatedState.Position( ).x ),
               convertKilometersToMeters( propagatedState.Position( ).y ),
               convertKilometersToMeters( propagatedState.Position( ).z ),
               convertKilometersToMeters( propagatedState.Velocity( ).x ),
               convertKilometersToMeters( propagatedState.Velocity( ).y ),
               convertKilometersToMeters( propagatedState.Velocity( ).z ) ).finished( );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cartesianState, recomputedCartesianState, 1.0e-8 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace atom
