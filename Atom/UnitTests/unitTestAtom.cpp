/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#define CATCH_CONFIG_MAIN

#include <catch.hpp>
 
namespace atom
{
namespace unit_tests
{

TEST_CASE( "Dummy", "[dummy]")
{
    REQUIRE( 1 == 1 );
}

} // namespace unit_tests
} // namespace atom
