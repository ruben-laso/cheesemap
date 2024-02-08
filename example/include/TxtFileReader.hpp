//
// Created by miguelyermo on 6/8/21.
// Modified by Ruben Laso on 7/2/24.
//

#pragma once

#include "FileReader.hpp"

#include <filesystem>
#include <fstream>

/**
 * @brief Specialization of FileRead to read .txt/.xyz files
 */
class TxtFileReader : public FileReader
{
	public:
	uint8_t numCols{};

	// ***  CONSTRUCTION / DESTRUCTION  *** //
	// ************************************ //
	TxtFileReader(const fs::path & path) : FileReader(path) {}
	~TxtFileReader() {}

	/**
	 * @brief Reads the points contained in the .txt/.xyz file
	 * @return Vector of Lpoint
	 */
	std::vector<chs::Point> read()
	{
		const auto split_line = [](std::string & line) {
			std::istringstream                 buf(line);
			std::istream_iterator<std::string> beg(buf), end;
			std::vector<std::string>           tokens(beg, end);

			return tokens;
		};

		std::ifstream file(path_.string());
		std::string   line{};

		if (not file.is_open()) { throw std::runtime_error("Could not open file"); }

		std::vector<chs::Point> points;

		while (std::getline(file, line, '\n'))
		{
			auto tokens = split_line(line);
			points.emplace_back(arma::vec3{ std::stod(tokens[0]),    // x
			                                std::stod(tokens[1]),    // y
			                                std::stod(tokens[2]) }); // z
		}

		std::cout << "Read points: " << points.size() << "\n";
		return points;
	}
};

std::vector<std::string> splitLine(std::string & line);