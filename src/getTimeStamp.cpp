#include "getTimeStamp.hpp"

std::string getTimeStamp()
{
	// Begin a timer for the output directory.
    auto start = std::chrono::system_clock::now();

	// Create string from the time the program started.
	time_t startTime = std::chrono::system_clock::to_time_t(start);
	std::string outputName = ctime(&startTime);

	// Strip out and replace difficult characters.
	std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ' ' ? '_' : ch;});
	std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ':' ? '-' : ch;});
	outputName.erase(std::remove(outputName.begin(), outputName.end(), '\n'), outputName.end());

	// Return the string so that the user can use it as appropriate.
	return outputName;

}