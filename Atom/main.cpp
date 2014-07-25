/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "Atom/convertCartesianStateToTwoLineElements.h"

//! Execute main-function.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
	const double datum = 1.234567;

	std::cout << atom::printElement( datum, 10, ' ' ) << std::endl;

	std::ofstream file("/home/kartikkumar/Applications/atom/bin/test.txt");

	file << atom::printElement( datum, 10, ' ' ) << std::endl;

	file.close( );

	return EXIT_SUCCESS;
}