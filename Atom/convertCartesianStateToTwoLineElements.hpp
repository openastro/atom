/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H
#define ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/Globals.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <SAM/sam.hpp>
#include <SML/sml.hpp>

#include <Atom/printFunctions.hpp>

namespace atom
{

//! Convert Cartesian state to TLE (Two Line Elements).
/*!
 * Converts a given Cartesian state (position, velocity) to an equivalent TLE. 
 *
 * This function makes use of a root-finder to solve a non-linear system. Locating the root of the 
 * non-linear system corresponds finding the TLE that, when evaluated at its epoch using the 
 * SGP4/SDP4 propagator (Vallado, 2012), yields the target Cartesian state (within tolerance).
 *
 * Details of the underlying non-linear system and algorithm are catalogued by 
 * Kumar, et al. (2014).
 *
 * @sa     evaluateCartesianToTwoLineElementsSystem, DateTime
 * @tparam Integer                     Type for integers
 * @tparam Real                        Type for reals
 * @tparam Vector                      Type for vector of reals
 * @param  cartesianState              Cartesian state [km; km/s]
 * @param  epoch                       Epoch associated with Cartesian state, stored in a
 *                                     DateTime object
 * @param  solverStatusSummary         Status of non-linear solver printed as a table
 * @param  numberOfIterations          Number of iterations completed by solver
 * @param  referenceTle                Reference Two Line Elements. This is used as reference to
 *                                     construct the effective TLE for the given Cartesian state
 *                                     [default: 0-TLE].
 * @param  earthGravitationalParameter Earth gravitational parameter [km^3 s^-2] [default: mu_SGP]
 * @param  absoluteTolerance           Absolute tolerance used to check if root-finder has
 *                                     converged [default: 1.0e-10] (see Kumar, et al. (2014) for
 *                                     details on how convergence is tested)
 * @param  relativeTolerance           Relative tolerance used to check if root-finder has
 *                                     converged [default: 1.0e-5] (see Kumar, et al. (2014) for
 *                                     details on how convergence is tested)
 * @param  maximumIterations           Maximum number of solver iterations permitted. Once the
 *                                     solver reaches this limit, the loop will be broken and the
 *                                     solver status will report that it has not converged
 *                                     [default: 100].
 * @return                             TLE object that generates target Cartesian state when
 *                                     propagated with SGP4 propagator to target epoch
 */
template< typename Integer, typename Real, typename Vector >
const Tle convertCartesianStateToTwoLineElements( 
    const Vector& cartesianState,
    const DateTime& epoch,
    std::string& solverStatusSummary,
    Integer& numberOfIterations,
    const Tle referenceTle = Tle( ),
    const Real earthGravitationalParameter = kMU,
    const Real absoluteTolerance = 1.0e-10,
    const Real relativeTolerance = 1.0e-5,
    const Integer maximumIterations = 100 );

//! Convert Cartesian state to TLE (Two Line Elements).
/*!
 * Converts a given Cartesian state (position, velocity) to an equivalent TLE. 
 *
 * This function makes use of a root-finder to solve a non-linear system. Locating the root of the 
 * non-linear system corresponds finding the TLE that, when evaluated at its epoch using the 
 * SGP4/SDP4 propagator (Vallado, 2012), yields the target Cartesian state (within tolerance).
 *
 * Details of the underlying non-linear system and algorithm are catalogued by 
 * Kumar, et al. (2014).
 *
 * This is a function overload to ensure that the user can opt to leave out solver summary status
 * string and number of iterations counter from the call (overload is necessary since non-const
 * references cannot be assigned default values in C++).
 *
 * @sa     convertCartesianStateToTwoLineElements, evaluateCartesianToTwoLineElementsSystem,
 *         DateTime
 * @tparam Integer                     Type for integers
 * @tparam Real                        Type for reals
 * @tparam Vector                      Type for vector of reals
 * @param  cartesianState              Cartesian state [km; km/s]
 * @param  epoch                       Epoch associated with Cartesian state, stored in a
 *                                     DateTime object
 * @return                             TLE object that generates target Cartesian state when
 *                                     propagated with SGP4 propagator to target epoch
 */
template< typename Integer, typename Real, typename Vector >
const Tle convertCartesianStateToTwoLineElements( 
    const Vector& cartesianState, const DateTime& epoch );

//! Compute residuals for converting Cartesian state to TLE.
/*!
 * Evaluates system of non-linear equations and computes residuals to find TLE corresponding with 
 * target Cartesian state. The residual function, \f$R\f$ is computed as follows:
 *  \f[ 
 *      R = 0 = \bar{x}^{new} - \bar{x}^{target}
 *  \f]
 * where \f$\bar{x}^{new}\f$ is the new Cartesian state computed by updating the TLE mean elements 
 * and propagating the TLE using the SGP4 propagator, and \f$\bar{x}^{target}\f$ is the target
 * Cartesian state. Note that the residuals are used to drive a root-finding process that uses the
 * GSL library.
 * 
 * @sa convertCartesianStateToTwoLineElements
 * @tparam Real                 Type for reals
 * @tparam Vector               Type for vector of reals
 * @param  independentVariables Vector of independent variables used by the root-finder
 * @param  parameters           Parameters required to compute the objective function
 * @param  residuals            Vector of computed residuals
 * @return                      GSL flag indicating success or failure
 */
template< typename Real, typename Vector >
int computeCartesianToTwoLineElementResiduals( const gsl_vector* independentVariables, 
                                               void* parameters, 
                                               gsl_vector* residuals );

//! Update TLE mean elements.
/*!
 * Updates mean elements stored in TLE based on current osculating elements. This function 
 * uses the osculating elements to replace mean elements (converts units and computes mean anomaly
 * and mean motion).
 *
 * @tparam Real                        Real type
 * @param  newKeplerianElements        New Keplerian elements
 * @param  oldTle                      TLE in which the mean elements are to be replaced
 * @param  earthGravitationalParameter Earth gravitational parameter [km^3 s^-2]
 * @return                             New TLE with mean elements updated.
 */
template< typename Real >
const Tle updateTleMeanElements( const gsl_vector* newKeplerianElements, 
                                 const Tle oldTle,
                                 const Real earthGravitationalParameter );

//! Parameter struct used by Cartesian-to-TLE residual function.
/*!
 * Data structure with parameters used to compute Cartesian-to-TLE residual function.
 *
 * @sa computeCartesianToTwoLineElementResiduals
 * @tparam Real   Type for reals
 * @tparam Vector Type for vector of reals
 */
template< typename Real, typename Vector >
struct CartesianToTwoLineElementsParameters;

//! Convert Cartesian state to TLE (Two Line Elements).
template< typename Integer, typename Real, typename Vector >
const Tle convertCartesianStateToTwoLineElements( 
    const Vector& cartesianState,
    const DateTime& epoch,
    std::string& solverStatusSummary,
    Integer& numberOfIterations,
    const Tle referenceTle,
    const Real earthGravitationalParameter,
    const Real absoluteTolerance,
    const Real relativeTolerance,
    const Integer maximumIterations )
{
    // Compute current state in Keplerian elements.
    Vector keplerianElements = sam::convertCartesianToKeplerianElements(
        cartesianState, earthGravitationalParameter );

    // Store reference TLE as new TLE and update epoch.
    Tle newTle = referenceTle;
    newTle.updateEpoch( epoch );  

    // Set up parameters for residual function.
    CartesianToTwoLineElementsParameters< Real, Vector > parameters( 
        cartesianState, earthGravitationalParameter, newTle ); 

    // Set up residual function.
    gsl_multiroot_function cartesianToTwoLineElementsFunction
        = { &computeCartesianToTwoLineElementResiduals< Real, Vector >, 
            6, 
            &parameters };

    // Set initial guess.
    gsl_vector* initialGuess = gsl_vector_alloc( 6 );
    for ( Integer i = 0; i < 6; i++ )
    {
        gsl_vector_set( initialGuess, i, keplerianElements[ i ] );      
    }

    // Set up solver type (derivative free).
    const gsl_multiroot_fsolver_type* solverType = gsl_multiroot_fsolver_hybrids;

    // Allocate memory for solver.
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc( solverType, 6 );

    // Set solver to use residual function with initial guess.
    gsl_multiroot_fsolver_set( solver, &cartesianToTwoLineElementsFunction, initialGuess );

     // Declare current solver status and iteration counter.
    Integer solverStatus = false;
    Integer counter = 0;

    // Set up buffer to store solver status summary table.
    std::ostringstream summary;

    // Print header for summary table to buffer.
    summary << printSolverStateTableHeader( );

    do
    {
        // Print current state of solver for summary table.
        summary << printSolverState( counter, solver );

        // Increment iteration counter.
        ++counter;

        // Execute solver iteration.
        solverStatus = gsl_multiroot_fsolver_iterate( solver );

        // Check if solver is stuck; if it is stuck, break from loop.
        if ( solverStatus )   
        {
            std::cerr << solverStatus << std::endl;
            std::cerr << summary.str( ) << std::endl;
            std::cerr << std::endl;
            throw std::runtime_error( "ERROR: Non-linear solver is stuck!" );
        }

        // Check if root has been found (within tolerance).
        solverStatus = gsl_multiroot_test_delta( 
          solver->dx, solver->x, absoluteTolerance, relativeTolerance );

    } while ( solverStatus == GSL_CONTINUE && counter < maximumIterations );

    // Save number of iterations.
    numberOfIterations = counter - 1;

    // Print final status of solver to buffer.
    summary << std::endl;
    summary << "Status of non-linear solver: " << gsl_strerror( solverStatus ) << std::endl;
    summary << std::endl;

    // Write buffer contents to solver status summary string.
    solverStatusSummary = summary.str( );

    // Update TLE with converged mean elements.
    newTle = updateTleMeanElements( solver->x, newTle, earthGravitationalParameter );

    // Free up memory.
    gsl_multiroot_fsolver_free( solver );
    gsl_vector_free( initialGuess );

    return newTle;
}

//! Convert Cartesian state to TLE (Two Line Elements).
template< typename Integer, typename Real, typename Vector >
const Tle convertCartesianStateToTwoLineElements( 
    const Vector& cartesianState, const DateTime& epoch )
{
    std::string dummyString = "";
    Integer dummyInteger = 0;
    return convertCartesianStateToTwoLineElements< Integer, Real >(
        cartesianState, epoch, dummyString, dummyInteger );
}

//! Compute residuals for converting Cartesian state to TLE.
template< typename Real, typename Vector >
int computeCartesianToTwoLineElementResiduals( const gsl_vector* independentVariables, 
                                               void* parameters, 
                                               gsl_vector* residuals )
{
    // Store reference TLE.
    const Tle referenceTle 
        = static_cast< CartesianToTwoLineElementsParameters< Real, Vector >* >( 
            parameters )->referenceTle;

    // Store gravitational parameter.
    const Real earthGravitationalParameter 
        = static_cast< CartesianToTwoLineElementsParameters< Real, Vector >* >( 
            parameters )->earthGravitationalParameter;

    // Update TLE mean elements and store as new TLE.
    Tle newTle = updateTleMeanElements( 
        independentVariables, referenceTle, earthGravitationalParameter );

    // Propagate new TLE to epoch of TLE.
    SGP4 sgp4( newTle );
    Eci propagatedState = sgp4.FindPosition( 0.0 );

    // Store target state.
    const Vector targetState
        = static_cast< CartesianToTwoLineElementsParameters< Real, Vector >* >( 
            parameters )->targetState;

    // Compute circular velocity at Earth radius (scaling fact used to non-dimensionalize 
    // velocities) [km/s].
    const Real circularVelocityEarthRadius = sam::computeCircularVelocity( kXKMPER, kMU );

    // Evaluate system of non-linear equations and store residuals.    
    gsl_vector_set( residuals, 0, ( propagatedState.Position( ).x - targetState[ 0 ] ) / kXKMPER );
    gsl_vector_set( residuals, 1, ( propagatedState.Position( ).y - targetState[ 1 ] ) / kXKMPER );
    gsl_vector_set( residuals, 2, ( propagatedState.Position( ).z - targetState[ 2 ] ) / kXKMPER );
    gsl_vector_set( residuals, 3, 
                    ( propagatedState.Velocity( ).x - targetState[ 3 ] ) 
                    / circularVelocityEarthRadius );
    gsl_vector_set( residuals, 4, 
                    ( propagatedState.Velocity( ).y - targetState[ 4 ] ) 
                    / circularVelocityEarthRadius );
    gsl_vector_set( residuals, 5, 
                    ( propagatedState.Velocity( ).z - targetState[ 5 ] ) 
                    / circularVelocityEarthRadius );

    return GSL_SUCCESS;        
}

//! Update TLE mean elements.
template< typename Real >
const Tle updateTleMeanElements( const gsl_vector* newKeplerianElements, 
                                 const Tle oldTle,
                                 const Real earthGravitationalParameter )
{
    // Copy old TLE to new object.
    Tle newTle( oldTle );
    
    // Compute new mean inclination [deg].
    const Real newInclination
        = sml::computeModulo( 
            sml::convertRadiansToDegrees( 
                gsl_vector_get( newKeplerianElements, sam::inclinationIndex ) ), 360.0 );

    // Compute new mean right ascending node [deg].
    const Real newRightAscendingNode 
        = sml::computeModulo( 
            sml::convertRadiansToDegrees( 
              gsl_vector_get( newKeplerianElements, 
                              sam::longitudeOfAscendingNodeIndex ) ), 360.0 );   

    // Compute new mean eccentricity [-].
    const Real newEccentricity 
        = gsl_vector_get( newKeplerianElements, sam::eccentricityIndex );

    // Compute new mean argument of perigee [deg].
    const Real newArgumentPerigee 
        = sml::computeModulo( 
            sml::convertRadiansToDegrees( 
                gsl_vector_get( newKeplerianElements, sam::argumentOfPeriapsisIndex ) ), 360.0 ); 

    // Compute new eccentric anomaly [rad].
    const Real eccentricAnomaly
        = sam::convertTrueAnomalyToEccentricAnomaly( 
            gsl_vector_get( newKeplerianElements, sam::trueAnomalyIndex ), newEccentricity );     

    // Compute new mean mean anomaly [deg].
    const Real newMeanAnomaly
        = sml::computeModulo( 
            sml::convertRadiansToDegrees( 
                sam::convertEccentricAnomalyToMeanAnomaly( 
                    eccentricAnomaly, newEccentricity ) ), 360.0 );

    // Compute new mean motion [rev/day].
    const Real newMeanMotion
        = sam::computeKeplerMeanMotion( 
            gsl_vector_get( newKeplerianElements, sam::semiMajorAxisIndex ),
            earthGravitationalParameter ) / ( 2.0 * sml::SML_PI ) * sam::SAM_JULIAN_DAY_IN_SECONDS;

    // Update mean elements in TLE with osculating elements.
    newTle.updateMeanElements( newInclination, 
                               newRightAscendingNode,
                               newEccentricity,
                               newArgumentPerigee,
                               newMeanAnomaly, 
                               newMeanMotion ); 

    return newTle;  
}

//! Parameter struct used by Cartesian-to-TLE residual function.
template< typename Real, typename Vector >
struct CartesianToTwoLineElementsParameters
{ 
public:

    //! Constructor taking parameter values.
    /*!
     * Default constructor, taking parameters for Cartesian-to-Two-Line-Elements conversion.
     * @sa convertCartesianStateToTwoLineElements, computeCartesianToTwoLineElementResiduals
     * @param aTargetState                  Target Cartesian state
     * @param anEarthGravitationalParameter Earth gravitational parameter [m^3 s^-2]
     * @param aReferenceTle                 Reference Two-Line-Elements
     */
    CartesianToTwoLineElementsParameters( 
        const Vector aTargetState,
        const Real anEarthGravitationalParameter,
        const Tle aReferenceTle )
        : targetState( aTargetState ),
          earthGravitationalParameter( anEarthGravitationalParameter ),
          referenceTle( aReferenceTle )
    { }

    //! Target state in Cartesian elements.
    const Vector targetState;

    //! Earth gravitational parameter [m^3 s^-2].
    const Real earthGravitationalParameter;

    //! Reference TLE.
    const Tle referenceTle;

protected:

private:
};

} // namespace atom

#endif // ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_H
