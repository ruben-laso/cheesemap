//
// Created by miguelyermo on 1/3/20.
//

/*
* FILENAME :  handlers.h  
* PROJECT  :  rule-based-classifier-cpp
* DESCRIPTION :
*  
*
*
*
*
* AUTHOR :    Miguel Yermo        START DATE : 03:07 1/3/20
*
*/

#ifndef CPP_HANDLERS_H
#define CPP_HANDLERS_H

#include <filesystem> // File extensions
#include <string>
#include <vector>

#include "TxtFileReader.hpp"

std::vector<chs::Point> readPointCloud(const std::filesystem::path & file)
{
	TxtFileReader fileReader{ file };

	return fileReader.read();
}

#endif //CPP_HANDLERS_H
