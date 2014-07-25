/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <stdexcept>

#include <libsgp4/Eci.h>  
#include <libsgp4/SGP4.h> 

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>  

#include "Atom/convertCartesianStateToTwoLineElements.h"

namespace atom
{

//! Convert Cartesian state to TLE (Two Line Elements).
const Tle convertCartesianStateToTwoLineElements( const Vector6d cartesianState,
                                                  const DateTime epoch,
                                                  const double earthGravitationalParameter,
                                                  const Tle referenceTle, 
                                                  std::string& solverStatusSummary,                                                 
                                                  const double tolerance,
                                                  const int maximumIterations )
{
    // Declare using-statements.
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;

    // Store reference TLE as new TLE.
    Tle newTle = referenceTle;

    // Set up SGP4 propagator.
    SGP4 sgp4( newTle );

    // Propagate new TLE to target epoch.
    Eci propagatedState = sgp4.FindPosition( epoch );

    // Store propagated state as current state.
    const Vector6d currentState
        = ( Eigen::VectorXd( 6 )
                << convertKilometersToMeters( propagatedState.Position( ).x ),
                   convertKilometersToMeters( propagatedState.Position( ).y ),
                   convertKilometersToMeters( propagatedState.Position( ).z ),
                   convertKilometersToMeters( propagatedState.Velocity( ).x ),
                   convertKilometersToMeters( propagatedState.Velocity( ).y ),
                   convertKilometersToMeters( propagatedState.Velocity( ).z ) ).finished( );

    // Compute current state in Keplerian elements.
    const Eigen::VectorXd currentStateInKeplerianElements
        = convertCartesianToKeplerianElements( currentState, earthGravitationalParameter );

    // Update TLE epoch.
    newTle.updateEpoch( epoch );

    // Set up parameters for objective function.
    CartesianToTwoLineElementsObjectiveParameters parameters( 
        cartesianState, earthGravitationalParameter, newTle );

    // Set up function.
    gsl_multiroot_function cartesianToTwoLineElementsFunction
        = {
             &evaluateCartesianToTwoLineElementsSystem, 
             6, 
             &parameters
          };

    // Set initial guess.
    gsl_vector* initialGuess = gsl_vector_alloc( 6 );
    for ( unsigned int i = 0; i < 6; i++ )
    {
        gsl_vector_set( initialGuess, i, currentStateInKeplerianElements[ i ] );      
    }

    // Set up solver type (derivative free).
    const gsl_multiroot_fsolver_type* solverType = gsl_multiroot_fsolver_hybrids;

    // Allocate memory for solver.
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc( solverType, 6 );

    // Set solver to use function with initial guess.
    gsl_multiroot_fsolver_set( solver, &cartesianToTwoLineElementsFunction, initialGuess );

    // Declare current solver status and iteration counter.
    int solverStatus;
    int iterationCounter = 0;

    // Set up buffer to store solver status summary table.
    std::ostringstream buffer;

    // Print header for summary table to buffer.
    buffer << printSolverStateTableHeader( );

    do
    {
        // Print current state of solver for summary table.
        buffer << printSolverState( iterationCounter, solver );        

        // Increment iteration counter.
        ++iterationCounter;

        // Execute solver iteration.
        solverStatus = gsl_multiroot_fsolver_iterate( solver );

        // Check if solver is stuck; if it is stuck, break from loop.
        if ( solverStatus )   
        {
            throw std::runtime_error( "ERROR: Non-linear solver is stuck!" );
        }

        // Check if root has been found (within tolerance).
        solverStatus = gsl_multiroot_test_residual( solver->f, tolerance );

    } while ( solverStatus == GSL_CONTINUE && iterationCounter < maximumIterations );

    // Print final status of solver to buffer.
    buffer << std::endl;
    buffer << "Status of non-linear solver: " << gsl_strerror( solverStatus ) << std::endl;
    buffer << std::endl;

    // Store final Keplerian elements.
    Vector6d finalKeplerianElements( 6 );
    for ( unsigned int i = 0; i < 6; i++ )
    {
        finalKeplerianElements( i ) = gsl_vector_get( solver->x, i );
    }

    // Free up memory.
    gsl_multiroot_fsolver_free( solver );
    gsl_vector_free( initialGuess );

    // Update and return new TLE object.
    return updateTleMeanElements( finalKeplerianElements, newTle, earthGravitationalParameter );
}

//! Evaluate system of non-linear equations for converting Cartesian state to TLE.
int evaluateCartesianToTwoLineElementsSystem( const gsl_vector* independentVariables, 
                                              void* parameters, 
                                              gsl_vector* functionValues )
{
    // // Declare using-statements.
    // using namespace tudat::basic_astrodynamics::unit_conversions;

    // // Store Keplerian elements in vector of independent variables.
    // Eigen::VectorXd currentStateInKeplerianElements( 6 );
    // for ( unsigned int i = 0; i < 6; i++ )
    // {
    //     currentStateInKeplerianElements[ i ] = gsl_vector_get( independentVariables, i );
    // }

    // // Store old TLE.
    // const Tle oldTle = static_cast< CartesianToTwoLineElementsObjectiveParameters* >( 
    //     parameters )->referenceTle;

    // // Store gravitational parameter.
    // const double earthGravitationalParameter 
    //     = static_cast< CartesianToTwoLineElementsObjectiveParameters* >( 
    //         parameters )->earthGravitationalParameter;
 
    // // Update TLE mean elements and store as new TLE.
    // Tle newTle = updateTleMeanElements( 
    //     currentStateInKeplerianElements, oldTle, earthGravitationalParameter );

    // // Propagate new TLE to epoch of TLE.
    // SGP4 sgp4( newTle );
    // Eci propagatedState = sgp4.FindPosition( 0.0 );

    // // Store new state.
    // const Eigen::VectorXd newState 
    //     = ( Eigen::VectorXd( 6 ) 
    //         << convertKilometersToMeters( propagatedState.Position( ).x ),
    //            convertKilometersToMeters( propagatedState.Position( ).y ),
    //            convertKilometersToMeters( propagatedState.Position( ).z ),
    //            convertKilometersToMeters( propagatedState.Velocity( ).x ),
    //            convertKilometersToMeters( propagatedState.Velocity( ).y ),
    //            convertKilometersToMeters( propagatedState.Velocity( ).z ) ).finished( );

    // // Store target state.
    // const Eigen::VectorXd targetState
    //     = static_cast< CartesianToTwoLineElementsObjectiveParameters* >( parameters )->targetState;

    // // Evaluate system of non-linear equations and store function values.
    // for ( unsigned int i = 0; i < 6; i++ )
    // {
    //     // gsl_vector_set( functionValues, i, ( newState[ i ] - targetState[ i ] ) 
    //     //                                    * ( newState[ i ] - targetState[ i ] ) );
    //     gsl_vector_set( functionValues, i, ( newState[ i ] - targetState[ i ] ) );        
    // }

    // Return success.
    return GSL_SUCCESS;
}

//! Print header for table containing summary of non-linear solver state.
std::string printSolverStateTableHeader( )
{
    std::ostringstream buffer;

    buffer << printElement( "#", 3, ' ' )
           << printElement( "a", 15, ' ' )
           << printElement( "e", 15, ' ' )
           << printElement( "i", 15, ' ' )
           << printElement( "AoP", 15, ' ' )
           << printElement( "RAAN", 15, ' ' )
           << printElement( "TA", 15, ' ' )
           << printElement( "f1", 15, ' ' )
           << printElement( "f2", 15, ' ' )
           << printElement( "f3", 15, ' ' )
           << printElement( "f4", 15, ' ' )
           << printElement( "f5", 15, ' ' )
           << printElement( "f6", 15, ' ' )
           << std::endl;

    return buffer.str( );
}

//! Print current state of non-linear solver for summary table.
std::string printSolverState( int iteration, gsl_multiroot_fsolver* solver )
{
    std::ostringstream buffer;

    buffer << printElement( iteration, 3, ' ' )
           << printElement( gsl_vector_get( solver->x, 0 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 1 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 2 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 3 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 4 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->x, 5 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 0 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 1 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 2 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 3 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 4 ), 15, ' ' )
           << printElement( gsl_vector_get( solver->f, 5 ), 15, ' ' )
           << std::endl;    

    return buffer.str( );
}

} // namespace atom
