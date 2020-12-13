#pragma once

#include <algorithm>
#include "Chromosome.h"
#include "Utils.h"

enum class MutationOperator
{
	saltation,
	inversion
};

class Mutation
{
public:
	static void SaltationMutation(Chromosome& chromosome)
	{
		const std::size_t chromosomeLength = chromosome.Length();
		const std::size_t firstPoint = GenerateRandomInt(0, chromosomeLength - 1);
		const std::size_t secondPoint = GenerateRandomInt(0, chromosomeLength - 1);

		std::size_t tmp = chromosome.GetGene(firstPoint);
		chromosome.SetGene(firstPoint, chromosome.GetGene(secondPoint));
		chromosome.SetGene(secondPoint, tmp);

		chromosome.CheckCorrectness();
	}

	static void InversionMutation(Chromosome& chromosome)
	{
		const std::size_t chromosomeLength = chromosome.Length();
		std::size_t firstPoint = GenerateRandomInt(0, chromosomeLength - 1);
		std::size_t secondPoint = GenerateRandomInt(0, chromosomeLength - 1);
		if (firstPoint > secondPoint)
		{
			std::swap(firstPoint, secondPoint);
		}
		const std::size_t interval = secondPoint - firstPoint;

		for (std::size_t i = 0; i < interval / 2; ++i)
		{
			std::size_t tmp = chromosome.GetGene(firstPoint + i);
			chromosome.SetGene(firstPoint + i, chromosome.GetGene(secondPoint - i));
			chromosome.SetGene(secondPoint - i, tmp);
		}

		chromosome.CheckCorrectness();
	}
};
