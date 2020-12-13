#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "Utils.h"

enum class InitialPopulationHeuristic
{
	random,
	nearestNeighbor
};

class Chromosome
{
public:
	using DistanceMatrix = std::vector<std::vector<double>>;

	Chromosome(
		const DistanceMatrix& distanceMatrix,
		const InitialPopulationHeuristic& heuristic = InitialPopulationHeuristic::random) :
		m_chromosomeLength(distanceMatrix.size()),
		m_distanceMatrix(distanceMatrix),
		m_fitness(-1.0),
		m_distance(0.0)
	{
		switch (heuristic)
		{
		case InitialPopulationHeuristic::random:
			{
				for (std::size_t i = 0; i < m_chromosomeLength; ++i)
				{
					m_chromosome.push_back(i);
				}
				std::shuffle(m_chromosome.begin(), m_chromosome.end(), std::mt19937{ std::random_device{}() });

				CalculateFitness();
			}
			break;
		case InitialPopulationHeuristic::nearestNeighbor:
			{
				std::vector<std::size_t> citySet;
				for (std::size_t i = 0; i < m_chromosomeLength; ++i)
				{
					citySet.push_back(i);
				}
				std::size_t firstCity = cityCount;
				m_chromosome.push_back(firstCity);				
				citySet.erase(std::remove(citySet.begin(), citySet.end(), firstCity), citySet.end());

				std::size_t secondCity = 0;
				double min = 0.0;
				while (!citySet.empty())
				{
					secondCity = citySet.front();
					min = m_distanceMatrix.at(firstCity).at(secondCity);
					for (std::size_t i = 1; i < citySet.size(); ++i)
					{
						if (distanceMatrix.at(firstCity).at(citySet.at(i)) < min)
						{
							min = distanceMatrix.at(firstCity).at(citySet.at(i));
							secondCity = citySet.at(i);
						}
					}
					m_chromosome.push_back(secondCity);
					citySet.erase(std::remove(citySet.begin(), citySet.end(), secondCity), citySet.end());
					firstCity = secondCity;
				}

				cityCount = (++cityCount) % m_chromosomeLength;
				CalculateFitness();
			}
			break;
		default:
			break;
		}
	}

	void SetGene(const std::size_t& pos, const std::size_t& gene)
	{
		m_chromosome.at(pos) = gene;
		CalculateFitness();
	}

	std::size_t GetGene(const std::size_t& pos) const
	{
		return m_chromosome.at(pos);
	}

	std::vector<std::size_t> GetChromosome() const
	{
		return m_chromosome;
	}

	double GetFitness() const
	{
		return m_fitness;
	}

	double GetDistance() const
	{
		return m_distance;
	}

	std::size_t Length() const
	{
		return m_chromosome.size();
	}

	bool ContainsGene(const std::size_t& gene) const
	{
		return std::find(m_chromosome.begin(), m_chromosome.end(), gene) != m_chromosome.end();
	}

	void CheckCorrectness()
	{
		if (std::unique(m_chromosome.begin(), m_chromosome.end()) != m_chromosome.end())
		{
			assert(false);
		}
	}

private:
	void CalculateFitness()
	{
		double distance = 0;
		for (std::size_t i = 0; i < m_chromosome.size() - 1; ++i)
		{
			distance += m_distanceMatrix.at(m_chromosome.at(i)).at(m_chromosome.at(i + 1));
		}
		distance += m_distanceMatrix.at(m_chromosome.back()).at(m_chromosome.front());
		m_fitness = 1 / distance;
		m_distance = distance;
	}

	static std::size_t cityCount;
	std::size_t m_chromosomeLength;
	DistanceMatrix m_distanceMatrix;
	std::vector<std::size_t> m_chromosome;
	double m_fitness;
	double m_distance;
};

std::size_t Chromosome::cityCount = 0;

std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome)
{
	std::string delim = "";
	for (auto& i : chromosome.GetChromosome())
	{
		os << delim << i + 1;
		delim = " ";
	}
	return os;
}

bool operator==(const Chromosome& lhs, const Chromosome& rhs)
{
	const std::size_t length = std::min(lhs.Length(), rhs.Length());
	for (std::size_t i = 0; i < length; ++i)
	{
		if (lhs.GetGene(i) != rhs.GetGene(i))
		{
			return false;
		}
	}
	return true;
}

bool operator!=(const Chromosome& lhs, const Chromosome& rhs)
{
	return !(lhs == rhs);
}

std::size_t HammingDistance(const Chromosome& a, const Chromosome& b)
{
	const std::size_t chromosomeLength = std::min(a.Length(), b.Length());
	std::size_t distance = 0;
	for (std::size_t i = 0; i < chromosomeLength; ++i)
	{
		if (a.GetGene(i) != b.GetGene(i))
		{
			++distance;
		}
	}
	return distance;
}
