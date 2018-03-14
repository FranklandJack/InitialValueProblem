#include "makeDirectory.hpp"

std::string makeDirectory(const std::string &name)
{

	// Create directory path from the string.
	boost::filesystem::path outPath = name;

	// If user calls program more than once a second so that directories will be overwritten append an index.
	for(int i = 2; boost::filesystem::exists(outPath) && i < 100; ++i)
	{
	    std::stringstream ss;
	    ss << outPath << i;
	    outPath = ss.str();
	}

	// Create the directory for output.
	boost::filesystem::create_directory(outPath);

	// Return the string so that the user can create files in the directory.
	return name;

}