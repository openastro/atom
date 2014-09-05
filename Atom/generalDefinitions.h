/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef ATOM_GENERAL_DEFINITIONS_H
#define ATOM_GENERAL_DEFINITIONS_H

namespace atom
{

//! Typedef for Vector6d.
typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//! Declare empty string.
static std::string emptyString = std::string( "" );    

} // namespace atom

#endif // ATOM_GENERAL_DEFINITIONS_H