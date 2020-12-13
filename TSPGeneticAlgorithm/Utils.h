#pragma once

#include <algorithm>
#include <random>
#include <string>

std::size_t GenerateRandomInt(const std::size_t& minValue, const std::size_t& maxValue)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> uid(minValue, maxValue);
	return uid(gen);
}

double GenerateRandomReal(const double& minValue, const double& maxValue)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> urd(minValue, maxValue);
	return urd(gen);
}

bool IsNumber(const std::string& number)
{
	return std::all_of(number.begin(), number.end(), ::isdigit);
}
