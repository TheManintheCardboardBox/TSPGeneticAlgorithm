#pragma once

#include <algorithm>
#include "Chromosome.h"
#include "Population.h"
#include "Utils.h"

enum class SelectionOperator
{
	fitnessProportionate,
	tournament
};

class Selection
{
public:
	static Chromosome FitnessProportionateSelection(Population& population)
	{
		population.Shuffle();
		double randomPosition = GenerateRandomReal(0.0, 1.0) * population.GetFitness();
		double total = 0;

		for (auto& i : population.GetPopulation())
		{
			total += i.GetFitness();
			if (total >= randomPosition)
			{
				return i;
			}
		}
		return population.GetPopulation().back();
	}

	static Chromosome TournamentSelection(Population& population, const std::size_t& tournamentSize)
	{
		population.Shuffle();
		std::vector<Chromosome> tournament;
		for (std::size_t i = 0; i < tournamentSize; ++i)
		{
			tournament.push_back(population.GetChromosome(GenerateRandomInt(0, population.Size() - 1)));
		}
		std::sort(tournament.begin(), tournament.end(), [](const Chromosome& a, const Chromosome& b)
		{
			return a.GetFitness() > b.GetFitness();
		});
		return tournament.front();
	}
};
