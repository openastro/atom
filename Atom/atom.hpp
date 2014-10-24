/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef ATOM_SOLVER_H
#define ATOM_SOLVER_H

#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h> 
#include <libsgp4/Tle.h>

#include "Atom/convertCartesianStateToTwoLineElements.hpp"

namespace atom
{

//! Execute Atom solver.
/*!
 * Executes Atom solver to find the transfer orbit connecting two positions. The epoch of the 
 * departure position and the Time-of-flight need to be specified.
 *
 * The Atom solver is an analog of the Lambert solver (Lancaster and Blanchard, 1969;
 * Gooding, 1990; Izzo, 2014), that aims to find the conic section that bridges two positions, at
 * given epochs, by using impulsive manveuvers (Delta-V maneuvers) at departure and arrival. The
 * Atom solver aims to solver a similar orbital transfer, subject to perturbations. The 
 * perturbations taken into account are those encoded in the SGP4/SDP4 propagators (Vallado, 2006).
 *
 * Since the Atom solver makes use fo the SGP4/SDP4 propagators, it can currently only solve for
 * perturbed transfers around the Earth. As a result, the Earth's gravitational parameter is fixed,
 * as specified by the SGP4/SDP4 propagators (Vallado, 2006).
 * 
 * Details of the underlying non-linear system and algorithm are catalogued by 
 * Kumar, et al. (2014).
 *
 * @sa     convertCartesianStateToTwoLineElements
 * @tparam Real                        Type for reals
 * @tparam Vector3                     Type for 3-vector of reals
 * @param  departurePosition           Cartesian position vector at departure [km]
 * @param  departureEpoch              Modified Julian Date (MJD) of departure
 * @param  arrivalPosition             Cartesian position vector at arrival [km]
 * @param  timeOfFlight                Time-of-Flight for orbital transfer [s]
 * @param  departureVelocityGuess      Initial guess for the departure velocity (serves as initial
 *                                     guess for the internal root-finding procedure) [km/s]
 * @param  solverStatusSummary         Status of non-linear solver printed as a table 
 * @param  numberOfIterations          Number of iterations completed by solver 
 * @param  referenceTle                Reference Two Line Elements [default: 0-TLE]
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
 * @return                             Departure and arrival velocities (stored in that order)
 */
template< typename Real, typename Vector3 >
const std::pair< Vector3, Vector3 > executeAtomSolver( 
    const Vector3& departurePosition, 
    const DateTime& departureEpoch,
    const Vector3& arrivalPosition, 
    const Real timeOfFlight, 
    const Vector3& departureVelocityGuess,
    std::string& solverStatusSummary,
    int& numberOfIterations,
    const Tle& referenceTle = Tle( ),
    const Real earthGravitationalParameter = kMU,
    const Real earthMeanRadius = kXKMPER,    
    const Real absoluteTolerance = 1.0e-10,
    const Real relativeTolerance = 1.0e-5,
    const int maximumIterations = 100 );

//! Execute Atom solver.
/*!
 * Executes Atom solver to find the transfer orbit connecting two positions. The epoch of the 
 * departure position and the Time-of-flight need to be specified.
 *
 * The Atom solver is an analog of the Lambert solver (Lancaster and Blanchard, 1969;
 * Gooding, 1990; Izzo, 2014), that aims to find the conic section that bridges two positions, at
 * given epochs, by using impulsive manveuvers (Delta-V maneuvers) at departure and arrival. The
 * Atom solver aims to solver a similar orbital transfer, subject to perturbations. The 
 * perturbations taken into account are those encoded in the SGP4/SDP4 propagators (Vallado, 2006).
 *
 * Since the Atom solver makes use fo the SGP4/SDP4 propagators, it can currently only solve for
 * perturbed transfers around the Earth. As a result, the Earth's gravitational parameter is fixed,
 * as specified by the SGP4/SDP4 propagators (Vallado, 2006).
 * 
 * Details of the underlying non-linear system and algorithm are catalogued by 
 * Kumar, et al. (2014).
 *
 * This is a function overload to ensure that the user can opt to leave out solver summary status
 * string and number of iterations counter from the call (overload is necessary since non-const
 * references cannot be assigned default values in C++).
 *
 * @sa     convertCartesianStateToTwoLineElements
 * @tparam Real                   Type for reals
 * @tparam Vector3                Type for 3-vector of reals
 * @param  departurePosition      Cartesian position vector at departure [km]
 * @param  departureEpoch         Modified Julian Date (MJD) of departure
 * @param  arrivalPosition        Cartesian position vector at arrival [km]
 * @param  timeOfFlight           Time-of-Flight for orbital transfer [s]
 * @param  departureVelocityGuess Initial guess for the departure velocity (serves as initial
 *                                guess for the internal root-finding procedure) [km/s]
 * @return                        Departure and arrival velocities (stored in that order)
 */
template< typename Real, typename Vector3 >
const std::pair< Vector3, Vector3 > executeAtomSolver( 
    const Vector3& departurePosition, 
    const DateTime& departureEpoch,
    const Vector3& arrivalPosition, 
    const Real timeOfFlight, 
    const Vector3& departureVelocityGuess );

//! Compute residuals to execute Atom solver.
/*!
 * Evaluates system of non-linear equations and computes residuals to execute the Atom solver. The
 * residual function, \f$\bar{R}\f$ is computed as follows:
 * The system of non-linear equations used is:
 *  \f[ 
 *      \bar{R} = 0 = \frac{\bar{r}_{new} - \bar{r}_{target}}{R_{Earth}}
 *  \f]
 * where \f$\bar{r}_{new}\f$ is the Cartesian position computed by propagating the 
 * initial, prescribed state under the action of an initial impulsive Delta V, by a prescribed 
 * time-of-flight, \f$\bar{r}_{target}\f$ is the target Cartesian position and 
 * \f$R_{Earth}\f$ is the mean radius of the Earth. Note that the residuals are non-dimensional.
 * They are used to drive a root-finding process that uses the GSL library.
 *
 * @sa executeAtomSolver
 * @tparam Real                 Type for reals
 * @param  independentVariables Vector of independent variables used by the root-finder
 * @param  parameters           Parameters required to compute the objective function
 * @param  residuals            Vector of computed residuals
 * @return                      GSL flag indicating success or failure
 */
template< typename Real >
int computeAtomResiduals( const gsl_vector* independentVariables,
                          void* parameters, 
                          gsl_vector* residuals );

//! Parameter struct used by Atom residual function.
/*!
 * Data structure with parameters used to compute Atom residual function.
 *
 * @sa computeAtomResiduals
 * @tparam Real    Type for reals
 * @tparam Vector3 Type for 3-vector of reals
 */
template< typename Real, typename Vector3 >
struct AtomParameters;

//! Execute Atom solver.
template< typename Real, typename Vector3 >
const std::pair< Vector3, Vector3 > executeAtomSolver( 
    const Vector3& departurePosition, 
    const DateTime& departureEpoch,
    const Vector3& arrivalPosition, 
    const Real timeOfFlight, 
    const Vector3& departureVelocityGuess,
    std::string& solverStatusSummary,
    int& numberOfIterations,
    const Tle& referenceTle,
    const Real earthGravitationalParameter,
    const Real earthMeanRadius,    
    const Real absoluteTolerance,
    const Real relativeTolerance,
    const int maximumIterations )
{
    // Set up parameters for residual function.
    AtomParameters< Real, Vector3 > parameters( departurePosition, 
                                                departureEpoch, 
                                                arrivalPosition, 
                                                timeOfFlight,
                                                earthGravitationalParameter, 
                                                earthMeanRadius, 
                                                referenceTle, 
                                                absoluteTolerance, 
                                                relativeTolerance, 
                                                maximumIterations ); 

    // Set up residual function.
    gsl_multiroot_function atomFunction = { &computeAtomResiduals< Real >, 3, &parameters };

    // Set initial guess.
    gsl_vector* initialGuess = gsl_vector_alloc( 3 );
    for ( int i = 0; i < 3; i++ )
    {
        gsl_vector_set( initialGuess, i, departureVelocityGuess[ i ] );      
    }

    // Set up solver type (derivative free).
    const gsl_multiroot_fsolver_type* solverType = gsl_multiroot_fsolver_hybrids;

    // Allocate memory for solver.
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc( solverType, 3 );

    // Set solver to use residual function with initial guess.
    gsl_multiroot_fsolver_set( solver, &atomFunction, initialGuess );

     // Declare current solver status and iteration counter.
    int solverStatus = false;
    int counter = 0;

    // Set up buffer to store solver status summary table.
    std::ostringstream summary;

    // Print header for summary table to buffer.
    summary << printAtomSolverStateTableHeader( );

    do
    {
        // Print current state of solver for summary table.
        summary << printAtomSolverState( counter, solver );

        // Increment iteration counter.
        ++counter;
        // Execute solver iteration.
        solverStatus = gsl_multiroot_fsolver_iterate( solver );

        // Check if solver is stuck; if it is stuck, break from loop.
        if ( solverStatus )   
        {
            std::cerr << "GSL solver status: " << solverStatus << std::endl;
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

    // Store final departure velocity.
    Vector3 departureVelocity( 3 );
    for ( int i = 0; i < 3; i++ )
    {
        departureVelocity[ i ] = gsl_vector_get( solver->x, i );
    }

    // Set departure state [km/s].
    std::vector< Real > departureState( 6 );
    for ( int i = 0; i < 3; i++ )
    {
        departureState[ i ] = departurePosition[ i ];
    }
    for ( int i = 0; i < 3; i++ )
    {
        departureState[ i + 3 ] = departureVelocity[ i ];
    }

    // Convert departure state to TLE.
    std::string dummyString = "";
    int dummyint = 0;
    const Tle departureTle = convertCartesianStateToTwoLineElements< Real >(
        departureState, 
        departureEpoch, 
        dummyString, 
        dummyint, 
        referenceTle, 
        earthGravitationalParameter, 
        earthMeanRadius, 
        absoluteTolerance, 
        relativeTolerance,
        maximumIterations );

    // Propagate departure TLE by time-of-flight using SGP4 propagator.
    SGP4 sgp4( departureTle );
    Eci arrivalState = sgp4.FindPosition( timeOfFlight );

    Vector3 arrivalVelocity( 3 );
    arrivalVelocity[ 0 ] = arrivalState.Velocity( ).x;
    arrivalVelocity[ 1 ] = arrivalState.Velocity( ).y;
    arrivalVelocity[ 2 ] = arrivalState.Velocity( ).z;

    // Free up memory.
    gsl_multiroot_fsolver_free( solver );
    gsl_vector_free( initialGuess );

    // Return departure and arrival velocities.
    return std::make_pair< Vector3, Vector3 >( departureVelocity, arrivalVelocity );
}

//! Execute Atom solver.
template< typename Real, typename Vector3 >
const std::pair< Vector3, Vector3 > executeAtomSolver( 
    const Vector3& departurePosition, 
    const DateTime& departureEpoch,
    const Vector3& arrivalPosition, 
    const Real timeOfFlight, 
    const Vector3& departureVelocityGuess )
{
    std::string dummyString = "";
    int dummyint = 0;
    return executeAtomSolver( departurePosition, 
                              departureEpoch, 
                              arrivalPosition, 
                              timeOfFlight, 
                              departureVelocityGuess,
                              dummyString, 
                              dummyint );
}

//! Compute residuals to execute Atom solver.
template< typename Real >
int computeAtomResiduals( const gsl_vector* independentVariables,
                          void* parameters, 
                          gsl_vector* residuals )
{
    // Store parameters locally.
    const std::vector< Real > departurePosition
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->departurePosition;

    const DateTime departureEpoch
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->departureEpoch;

    const std::vector< Real > targetPosition
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->targetPosition;

    const Real timeOfFlight
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->timeOfFlight;

    const Real earthGravitationalParameter
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->earthGravitationalParameter;

    const Real earthMeanRadius
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->earthMeanRadius;

    const Tle referenceTle
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->referenceTle;

    const Real absoluteTolerance
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->absoluteTolerance;

    const Real relativeTolerance
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->relativeTolerance;

    const int maximumIterations
        = static_cast< AtomParameters< Real, std::vector< Real > >* >( 
            parameters )->maximumIterations;

    // Set Departure state [km; km/s].
    std::vector< Real > departureVelocity( 3 );
    for ( int i = 0; i < 3; i++ )
    {
        departureVelocity[ i ] = gsl_vector_get( independentVariables, i );
    }

    std::vector< Real > departureState( 6 );
    for ( int i = 0; i < 3; i++ )
    {
        departureState[ i ] = departurePosition[ i ];
    }
    for ( int i = 0; i < 3; i++ )
    {
        departureState[ i + 3 ] = departureVelocity[ i ];
    }

    // Convert departure state to TLE.
    std::string dummyString = "";
    int dummyint = 0;
    const Tle departureTle = convertCartesianStateToTwoLineElements(
        departureState, 
        departureEpoch, 
        dummyString, 
        dummyint, 
        referenceTle, 
        earthGravitationalParameter, 
        earthMeanRadius, 
        absoluteTolerance, 
        relativeTolerance,
        maximumIterations );

    // Propagate departure TLE by time-of-flight using SGP4 propagator.
    SGP4 sgp4( departureTle );
    Eci arrivalState = sgp4.FindPosition( timeOfFlight );

    // Evaluate system of non-linear equations and store residuals.    
    gsl_vector_set( residuals, 0,
                    ( arrivalState.Position( ).x - targetPosition[ 0 ] ) / earthMeanRadius );
    gsl_vector_set( residuals, 1, 
                    ( arrivalState.Position( ).y - targetPosition[ 1 ] ) / earthMeanRadius );
    gsl_vector_set( residuals, 2, 
                    ( arrivalState.Position( ).z - targetPosition[ 2 ] ) / earthMeanRadius );

    return GSL_SUCCESS;        
}

//! Parameter struct used by Atom residual function.
template< typename Real, typename Vector3 >
struct AtomParameters
{ 
public:

    //! Constructor taking parameter values.
    /*!
     * Default constructor, taking parameters to execute Atom solver.
     * @sa executeAtomSolver, computeCartesianToTwoLineElementResiduals
     * @param aDeparturePosition            Cartesian departure position [km]
     * @param aDepartureEpoch               Modified Julian Date (MJD) of departure
     * @param aTargetPosition               Target Cartesian position [km]
     * @param aTimeOfFlight                 Time-of-Flight (TOF) [s]
     * @param anEarthGravitationalParameter Earth gravitational parameter [km^3 s^-2]
     * @param anEarthMeanRadius             Earth mean radius [km]     
     * @param aReferenceTle                 Reference Two-Line-Elements
     * @param anAbsoluteTolerance           Absolute tolerance used to check for convergence
     * @param aRelativeTolerance            Relative tolerance used to check for convergence
     * @param someMaximumIterations         Maximum number of solver iterations permitted
     */
    AtomParameters(
        const Vector3& aDeparturePosition,
        const DateTime& aDepartureEpoch,
        const Vector3& aTargetPosition,
        const Real aTimeOfFlight,
        const Real anEarthGravitationalParameter,
        const Real anEarthMeanRadius,        
        const Tle& aReferenceTle,
        const Real anAbsoluteTolerance,
        const Real aRelativeTolerance,
        const int someMaximumIterations )
        : departurePosition( aDeparturePosition ),
          departureEpoch( aDepartureEpoch ),
          targetPosition( aTargetPosition ),
          timeOfFlight( aTimeOfFlight ),
          earthGravitationalParameter( anEarthGravitationalParameter ),
          earthMeanRadius( anEarthMeanRadius ),          
          referenceTle( aReferenceTle ),
          absoluteTolerance( anAbsoluteTolerance ),
          relativeTolerance( aRelativeTolerance ),
          maximumIterations( someMaximumIterations )
    { }

    //! Departure position in Cartesian elements [km].
    const Vector3 departurePosition;

    //! Departure epoch in Modified Julian Date (MJD).
    const DateTime departureEpoch;

    //! Target position in Cartesian elements [km].
    const Vector3 targetPosition;

    //! Time-of-Flight (TOF) [s].
    const Real timeOfFlight;

    //! Earth gravitational parameter [km^3 s^-2].
    const Real earthGravitationalParameter;

    //! Earth mean radius [km].
    const Real earthMeanRadius;    

    //! Reference TLE.
    const Tle referenceTle;

    //! Absolute tolerance [-].
    const Real absoluteTolerance;

    //! Relative tolerance [-].
    const Real relativeTolerance;

    //! Maximum number of iterations.
    const int maximumIterations;

protected:

private:
};

} // namespace atom

#endif // ATOM_SOLVER_H
