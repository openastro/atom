/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_HPP
#define ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_HPP

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/Globals.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <Astro/astro.hpp>
#include <SML/sml.hpp>

#include <Atom/printFunctions.hpp>

namespace atom
{

const static int meanMotionIndex  = astro::semiMajorAxisIndex;
const static int meanAnomalyIndex = astro::trueAnomalyIndex;

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
 * @tparam Real                        Type for reals
 * @tparam Vector6                     Type for 6-vector of reals
 * @param  cartesianState              Cartesian state [km; km/s]
 * @param  epoch                       Epoch associated with Cartesian state, stored in a
 *                                     DateTime object
 * @param  solverStatusSummary         Status of non-linear solver printed as a table
 * @param  numberOfIterations          Number of iterations completed by solver
 * @param  referenceTle                Reference Two Line Elements. This is used as reference to
 *                                     construct the effective TLE for the given Cartesian state
 *                                     [default: 0-TLE].
 * @param  earthGravitationalParameter Earth gravitational parameter [km^3 s^-2] [default: mu_SGP]
 * @param  earthMeanRadius             Earth mean radius [km] [default: R_SGP]
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
template< typename Real, typename Vector6 >
const Tle convertCartesianStateToTwoLineElements(
    const Vector6& cartesianState,
    const DateTime& epoch,
    std::string& solverStatusSummary,
    int& numberOfIterations,
    const Tle& referenceTle = Tle( ),
    const Real earthGravitationalParameter = kMU,
    const Real earthMeanRadius = kXKMPER,
    const Real absoluteTolerance = 1.0e-10,
    const Real relativeTolerance = 1.0e-8,
    const int maximumIterations = 100 );

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
 * @tparam Real                        Type for reals
 * @tparam Vector6                     Type for 6-vector of reals
 * @param  cartesianState              Cartesian state [km; km/s]
 * @param  epoch                       Epoch associated with Cartesian state, stored in a
 *                                     DateTime object
 * @return                             TLE object that generates target Cartesian state when
 *                                     propagated with SGP4 propagator to target epoch
 */
template< typename Real, typename Vector6 >
const Tle convertCartesianStateToTwoLineElements( const Vector6& cartesianState,
                                                  const DateTime& epoch );

//! Compute residuals for converting Cartesian state to TLE.
/*!
 * Evaluates system of non-linear equations and computes residuals to find TLE corresponding with
 * target Cartesian state. The residual function, \f$\bar{R}\f$ is computed as follows:
 *  \f[
 *      \bar{R} = 0 = \begin{pmatrix}
 *                      \frac{\bar{r}_{new} - \bar{r}_{target}}{R_{Earth}}\\
 *                      \frac{\bar{v}_{new} - \bar{v}_{target}}{V_{c,Earth}}\\
 *                    \end{pmatrix}
 *  \f]
 * where \f$\bar{r}_{new}\f$ and \f$\bar{v}_{new}\f$ are the new Cartesian position and velocity
 * vectors, computed by updating the TLE mean elements and propagating the TLE using the SGP4
 * propagator, \f$\bar{r}_{target}\f$ anf \f$\bar{v}_{target}\f$ are the target Cartesian position
 * and velocity vectors, \f$R_{Earth}\f$ is the mean radius of the Earth, and \f$V_{c,Earth}\f$ is
 * the circular velocity at \f$R_{Earth}\f$. Note that the residuals are non-dimensional. They are
 * used to drive a root-finding process that uses the GSL library.
 *
 * @sa convertCartesianStateToTwoLineElements
 * @tparam Real                 Type for reals
 * @tparam Vector6              Type for 6-vector of reals
 * @param  independentVariables Vector of independent variables used by the root-finder
 * @param  parameters           Parameters required to compute the objective function
 * @param  residuals            Vector of computed residuals
 * @return                      GSL flag indicating success or failure
 */
template< typename Real, typename Vector6 >
int computeCartesianToTwoLineElementResiduals( const gsl_vector* independentVariables,
                                               void* parameters,
                                               gsl_vector* residuals );

//! Compute initial guess for TLE mean elements.
/*
 * Computes initial guess for TLE mean elements by converting Keplerian elements provided by user to
 * TLE mean elements. The TLE mean elements generated are to floating-point precision, since the
 * TLE class stores the internal values as doubles.
 *
 * @sa Tle
 * @tparam Real                         Type for reals
 * @tparam Vector6                      Type for 6-vector of reals
 * @param  keplerianElements            Keplerian elements to be used as initial guess
 *                                      The order of Keplerian elements in the vector is:
 *                                      - semi-major axis                                      [km]
 *                                      - eccentricity                                          [-]
 *                                      - inclination                                         [rad]
 *                                      - argument of periapsis                               [rad]
 *                                      - right ascension of ascending node                   [rad]
 *                                      - true anomaly                                        [rad]
 * @param  earthGravitationalParameter  Earth gravitational parameter                   [km^3 s^-2]
 * @return                              6-vector containing initial guess for TLE mean
 *                                      elements
 *                                      The order of mean TLE elements in the vector is:
 *                                      - mean inclination                                    [deg]
 *                                      - mean right ascension of ascending node              [deg]
 *                                      - mean eccentricity                                     [-]
 *                                      - mean argument of perigee                            [deg]
 *                                      - mean mean anomaly                                   [deg]
 *                                      - mean mean motion                                [rev/day]
 */
template< typename Real, typename Vector6 >
const Vector6 computeInitialGuessTleMeanElements( const Vector6& keplerianElements,
                                                  const Real earthGravitationalParameter );

//! Parameter struct used by Cartesian-to-TLE residual function.
/*!
 * Data structure with parameters used to compute Cartesian-to-TLE residual function.
 *
 * @sa computeCartesianToTwoLineElementResiduals
 * @tparam Vector6 Type for 6-vector of reals
 */
template< typename Vector6 >
struct CartesianToTwoLineElementsParameters;

//! Convert Cartesian state to TLE (Two Line Elements).
template< typename Real, typename Vector6 >
const Tle convertCartesianStateToTwoLineElements(
    const Vector6& cartesianState,
    const DateTime& epoch,
    std::string& solverStatusSummary,
    int& numberOfIterations,
    const Tle& referenceTle,
    const Real earthGravitationalParameter,
    const Real earthMeanRadius,
    const Real absoluteTolerance,
    const Real relativeTolerance,
    const int maximumIterations )
{
    // Store reference TLE as the template TLE and update epoch.
    Tle templateTle = referenceTle;
    templateTle.updateEpoch( epoch );

    // Set up parameters for residual function.
    CartesianToTwoLineElementsParameters< Vector6 > parameters( cartesianState, templateTle );

    // Set up residual function.
    gsl_multiroot_function cartesianToTwoLineElementsFunction
        = { &computeCartesianToTwoLineElementResiduals< Real, Vector6 >,
            6,
            &parameters };

    // Compute current state in Keplerian elements, for use as initial guess for the TLE mean
    // elements.
    const Vector6 initialKeplerianElements = astro::convertCartesianToKeplerianElements(
        parameters.targetState, earthGravitationalParameter );

    // Compute initial guess for TLE mean elements.
    const Vector6 initialTleMeanElements
        = computeInitialGuessTleMeanElements( initialKeplerianElements,
                                              earthGravitationalParameter );

    // Set initial guess.
    gsl_vector* initialGuessTleMeanElements = gsl_vector_alloc( 6 );
    for ( int i = 0; i < 6; i++ )
    {
        gsl_vector_set( initialGuessTleMeanElements, i, initialTleMeanElements[ i ] );
    }

    // Set up solver type (derivative free).
    const gsl_multiroot_fsolver_type* solverType = gsl_multiroot_fsolver_hybrids;

    // Allocate memory for solver.
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc( solverType, 6 );

    // Set solver to use residual function with initial guess for TLE mean elements.
    gsl_multiroot_fsolver_set( solver,
                               &cartesianToTwoLineElementsFunction,
                               initialGuessTleMeanElements );

     // Declare current solver status and iteration counter.
    int solverStatus = false;
    int counter = 0;

    // Set up buffer to store solver status summary table.
    std::ostringstream summary;

    // Print header for summary table to buffer.
    summary << printCartesianToTleSolverStateTableHeader( );

    do
    {
        // Print current state of solver for summary table.
        summary << printCartesianToTleSolverState( counter, solver );

        // Increment iteration counter.
        ++counter;

        // Execute solver iteration.
        solverStatus = gsl_multiroot_fsolver_iterate( solver );

        // Check if solver is stuck; if it is stuck, break from loop.
        if ( solverStatus )
        {
            solverStatusSummary = summary.str( );
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

    // Generate TLE with converged mean elements.
    Tle virtualTle = templateTle;

    Real convergedMeanEccentricity = gsl_vector_get( solver->x, 2 );
    if ( convergedMeanEccentricity < 0.0 )
    {
        convergedMeanEccentricity = std::fabs( gsl_vector_get( solver->x, 2 ) );
    }

    if ( convergedMeanEccentricity > 0.999 )
    {
        convergedMeanEccentricity = 0.99;
    }

    virtualTle.updateMeanElements( sml::computeModulo( std::fabs( gsl_vector_get( solver->x, 0 ) ), 180.0 ),
                                   sml::computeModulo( gsl_vector_get( solver->x, 1 ), 360.0 ),
                                   convergedMeanEccentricity,
                                   sml::computeModulo( gsl_vector_get( solver->x, 3 ), 360.0 ),
                                   sml::computeModulo( gsl_vector_get( solver->x, 4 ), 360.0 ),
                                   std::fabs( gsl_vector_get( solver->x, 5 ) ) );

    // Free up memory.
    gsl_multiroot_fsolver_free( solver );
    gsl_vector_free( initialGuessTleMeanElements );

    return virtualTle;
}

//! Convert Cartesian state to TLE (Two Line Elements).
template< typename Real, typename Vector6 >
const Tle convertCartesianStateToTwoLineElements( const Vector6& cartesianState,
                                                  const DateTime& epoch )
{
    std::string dummyString = "";
    int dummyint = 0;
    return convertCartesianStateToTwoLineElements< Real >(
        cartesianState, epoch, dummyString, dummyint );
}

//! Compute residuals for converting Cartesian state to TLE.
template< typename Real, typename Vector6 >
int computeCartesianToTwoLineElementResiduals( const gsl_vector* independentVariables,
                                               void* parameters,
                                               gsl_vector* residuals )
{
    const Vector6 targetState
        = static_cast< CartesianToTwoLineElementsParameters< Vector6 >* >( parameters )->targetState;
    const Tle templateTle
        = static_cast< CartesianToTwoLineElementsParameters< Vector6 >* >( parameters )->templateTle;

    // Create a TLE object with the mean TLE elements generated by the root-finder.
    Tle tle = templateTle;

    Real meanInclination = sml::computeModulo( std::fabs( gsl_vector_get( independentVariables, 0 ) ), 180.0 );
    Real meanRightAscendingNode = sml::computeModulo( gsl_vector_get( independentVariables, 1 ), 360.0 );

    Real meanEccentricity = gsl_vector_get( independentVariables, 2 );
    if ( meanEccentricity < 0.0 )
    {
        meanEccentricity = std::fabs( gsl_vector_get( independentVariables, 2 ) );
    }

    if ( meanEccentricity > 0.999 )
    {
        meanEccentricity = 0.99;
    }

    Real meanArgumentPerigee =  sml::computeModulo( gsl_vector_get( independentVariables, 3 ), 360.0 );
    Real meanMeanAnomaly =  sml::computeModulo( gsl_vector_get( independentVariables, 4 ), 360.0 );
    Real meanMeanMotion = std::fabs( gsl_vector_get( independentVariables, 5 ) );

    tle.updateMeanElements( meanInclination,
                            meanRightAscendingNode,
                            meanEccentricity,
                            meanArgumentPerigee,
                            meanMeanAnomaly,
                            meanMeanMotion );

    // Propagate the TLE object to the specified epoch using the SGP4 propagator.
    SGP4 sgp4( tle );
    Eci cartesianState = sgp4.FindPosition( 0.0 );

    // Compute residuals by computing the difference between the Cartesian state generated by the
    // SGP4 propagator and the target Cartesian state.
    gsl_vector_set( residuals, astro::xPositionIndex, cartesianState.Position( ).x - targetState[ astro::xPositionIndex ] );
    gsl_vector_set( residuals, astro::yPositionIndex, cartesianState.Position( ).y - targetState[ astro::yPositionIndex ] );
    gsl_vector_set( residuals, astro::zPositionIndex, cartesianState.Position( ).z - targetState[ astro::zPositionIndex ] );
    gsl_vector_set( residuals, astro::xVelocityIndex, cartesianState.Velocity( ).x - targetState[ astro::xVelocityIndex ] );
    gsl_vector_set( residuals, astro::yVelocityIndex, cartesianState.Velocity( ).y - targetState[ astro::yVelocityIndex ] );
    gsl_vector_set( residuals, astro::zVelocityIndex, cartesianState.Velocity( ).z - targetState[ astro::zVelocityIndex ] );

    return GSL_SUCCESS;
}

//! Compute initial guess for TLE mean elements.
template< typename Real, typename Vector6 >
const Vector6 computeInitialGuessTleMeanElements( const Vector6& keplerianElements,
                                                  const Real earthGravitationalParameter )

{
    Vector6 tleMeanElements = keplerianElements;

    // Compute mean inclination [deg].
    tleMeanElements[ 0 ]
    = sml::computeModulo(
        sml::convertRadiansToDegrees( keplerianElements[ astro::inclinationIndex ] ), 180.0 );

    // Compute mean right ascending node [deg].
    tleMeanElements[ 1 ]
        = sml::computeModulo(
            sml::convertRadiansToDegrees(
                keplerianElements[ astro::longitudeOfAscendingNodeIndex ] ), 360.0 );

    // Compute mean eccentricity [-].
    tleMeanElements[ 2 ] = keplerianElements[ astro::eccentricityIndex ];

    // Compute mean argument of perigee [deg].
    tleMeanElements[ 3 ]
        = sml::computeModulo(
            sml::convertRadiansToDegrees(
                keplerianElements[ astro::argumentOfPeriapsisIndex ] ), 360.0 );

    // Compute mean eccentric anomaly [rad].
    const Real eccentricAnomaly
        = astro::convertTrueAnomalyToEccentricAnomaly(
            keplerianElements[ astro::trueAnomalyIndex ],
            keplerianElements[ astro::eccentricityIndex ] );

    // Compute mean mean anomaly [deg].
    tleMeanElements[ 4 ]
        = sml::computeModulo(
            sml::convertRadiansToDegrees(
                astro::convertEccentricAnomalyToMeanAnomaly(
                    eccentricAnomaly, keplerianElements[ astro::eccentricityIndex ] ) ), 360.0 );

    // Compute new mean motion [rev/day].
    tleMeanElements[ 5 ]
        = astro::computeKeplerMeanMotion(
            keplerianElements[ astro::semiMajorAxisIndex ],
            earthGravitationalParameter ) / ( 2.0 * sml::SML_PI ) * astro::ASTRO_JULIAN_DAY_IN_SECONDS;

    return tleMeanElements;
}

//! Parameter struct used by Cartesian-to-TLE residual function.
template< typename Vector6 >
struct CartesianToTwoLineElementsParameters
{
public:

    //! Constructor taking parameter values.
    /*!
     * Default constructor, taking parameters for Cartesian-to-Two-Line-Elements conversion.
     * @sa convertCartesianStateToTwoLineElements, computeCartesianToTwoLineElementResiduals
     *
     * @tparam Vector6                      Vector of length 6
     * @param  aTargetState                 Target Cartesian state [km; km/s]
     * @param  aTemplateTle                 Template for Two-Line-Elements (TLE) with pre-filled
     *                                      values
     */
    CartesianToTwoLineElementsParameters(
        const Vector6& aTargetState,
        const Tle& aTemplateTle )
        : targetState( aTargetState ),
          templateTle( aTemplateTle )
    { }

    //! Target state in Cartesian elements [km; km/s].
    const Vector6 targetState;

    //! Template TLE that contains pre-filled values (e.g., epoch, Bstar).
    const Tle templateTle;

protected:

private:
};

} // namespace atom

#endif // ATOM_CONVERT_CARTESIAN_STATE_TO_TWO_LINE_ELEMENTS_HPP
