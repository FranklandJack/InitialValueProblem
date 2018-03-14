#ifndef makeDirectory_hpp
#define makeDirectory_hpp

#include <boost/filesystem.hpp> 
#include <string>
#include <sstream>
#include "getTimeStamp.hpp"

/**
 *\file
 *\brief function to create a directory.
 *\param constant string reference that is the name of the output directory.
 *\return string representing the name of the directory.
 *
 * The function will append integers to directory name if it already exists.
 */
std::string makeDirectory(const std::string &name);

#endif /* makeDirectory_hpp */

