#pragma once

#include <vector>
#include "Chromosome.h"
#include "Population.h"

enum class BreedingSystem
{
	panmixia,
	inbreeding,
	autobreeding
};

class Breeding
{
public:
	static std::vector<Chromosome> Panmixia(const Population& population)
	{
		const std::size_t populationSize = population.Size();
		const Chromosome parent1 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));
		Chromosome parent2 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));

		while (HammingDistance(parent1, parent2) == 0)
		{
			parent2 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));
		}

		return std::vector<Chromosome>{ parent1, parent2 };
	}

	static std::vector<Chromosome> Inbreeding(const Population& population, const std::size_t& inbreedingValue)
	{
		const std::size_t populationSize = population.Size();
		const Chromosome parent1 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));
		Chromosome parent2 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));

		while (HammingDistance(parent1, parent2) > inbreedingValue)
		{
			parent2 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));
		}

		return std::vector<Chromosome>{ parent1, parent2 };
	}

	static std::vector<Chromosome> Autobreeding(const Population& population, const std::size_t& autobreedingValue)
	{
		const std::size_t populationSize = population.Size();
		const Chromosome parent1 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));
		Chromosome parent2 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));

		while (HammingDistance(parent1, parent2) < autobreedingValue)
		{
			parent2 = population.GetChromosome(GenerateRandomInt(0, populationSize - 1));
		}

		return std::vector<Chromosome>{ parent1, parent2 };
	}
};
