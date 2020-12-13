#pragma once

#include <algorithm>
#include <vector>
#include "Chromosome.h"

class Population
{
public:
	using DistanceMatrix = std::vector<std::vector<double>>;

	Population() :
		m_populationSize(0),
		m_distanceMatrix(),
		m_population(),
		m_fitness(0)
	{
	}

	Population(
		const std::size_t& populationSize,
		const DistanceMatrix& distanceMatrix,
		const InitialPopulationHeuristic& heuristic = InitialPopulationHeuristic::random) :
		m_populationSize(populationSize),
		m_distanceMatrix(distanceMatrix),
		m_fitness(0)
	{
		switch (heuristic)
		{
		case InitialPopulationHeuristic::random:
			{
				while (m_population.size() < m_populationSize)
				{
					const Chromosome candidate = Chromosome(distanceMatrix, heuristic);
					if (ContainsChromosome(candidate))
					{
						m_population.push_back(Chromosome(distanceMatrix, heuristic));
					}					
				}
			}
			break;
		case InitialPopulationHeuristic::nearestNeighbor:
			{
				while (m_population.size() < m_populationSize)
				{
					m_population.push_back(Chromosome(distanceMatrix, heuristic));
				}
			}
			break;
		default:
			break;
		}
		std::shuffle(m_population.begin(), m_population.end(), std::mt19937{ std::random_device{}() });
		CalculateFitness();
	}

	void Add(const Chromosome& chromosome)
	{
		m_population.push_back(chromosome);
		m_populationSize++;
		CalculateFitness();
	}

	void OrderByFitness()
	{
		std::sort(m_population.begin(), m_population.end(), [](const Chromosome& a, const Chromosome& b)
		{
			return a.GetFitness() < b.GetFitness();
		});
	}

	void Shuffle()
	{
		std::shuffle(m_population.begin(), m_population.end(), std::mt19937{ std::random_device{}() });
	}

	Chromosome GetChromosome(const std::size_t& pos) const
	{
		return m_population.at(pos);
	}

	std::vector<Chromosome> GetPopulation() const
	{
		return m_population;
	}

	Chromosome GetFittest() const
	{
		std::vector<Chromosome> tmp = m_population;
		std::sort(tmp.begin(), tmp.end(), [](const Chromosome& a, const Chromosome& b)
		{
			return a.GetFitness() > b.GetFitness();
		});
		return tmp.front();
	}

	double GetFitness() const
	{
		return m_fitness;
	}

	std::size_t Size() const
	{
		return m_population.size();
	}

	bool ContainsChromosome(const Chromosome& chromosome)
	{
		return std::find(m_population.begin(), m_population.end(), chromosome) == m_population.end();
	}

private:
	void CalculateFitness()
	{
		double fitness = 0;
		for (auto& i : m_population)
		{
			fitness += i.GetFitness();
		}
		m_fitness = fitness;
	}

	std::size_t m_populationSize;
	DistanceMatrix m_distanceMatrix;
	std::vector<Chromosome> m_population;
	double m_fitness;
};
