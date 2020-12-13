#pragma once

#include <algorithm>
#include <iterator>
#include "Chromosome.h"
#include "Utils.h"

enum class CrossoverOperator
{
	OX,
	CX
};

class Crossover
{
public:
	static std::vector<Chromosome> OrderCrossover(const Chromosome& parent1, const Chromosome& parent2)
	{
		const std::size_t chromosomeLength = std::min(parent1.Length(), parent2.Length());
		Chromosome child1 = parent1;
		Chromosome child2 = parent2;

		const std::size_t firstPoint = GenerateRandomInt(1, chromosomeLength - 1);
		const std::size_t secondPoint = GenerateRandomInt(firstPoint, chromosomeLength - 1);

		for (std::size_t i = firstPoint; i < secondPoint; ++i)
		{
			child1.SetGene(i, parent1.GetGene(i));
			child2.SetGene(i, parent2.GetGene(i));
		}

		for (std::size_t i = secondPoint; ; i = (i + 1) % chromosomeLength)
		{
			for (std::size_t j = i; ; j = (++j) % chromosomeLength)
			{
				if (!child1.ContainsGene(parent2.GetGene(j)))
				{
					child1.SetGene(i, parent2.GetGene(j));
					break;
				}
				if (j == secondPoint - 1)
				{
					break;
				}
			}
			if (i == firstPoint - 1)
			{
				break;
			}
		}

		for (std::size_t i = secondPoint; ; i = (i + 1) % chromosomeLength)
		{
			for (std::size_t j = i; ; j = (++j) % chromosomeLength)
			{
				if (!child2.ContainsGene(parent1.GetGene(j)))
				{
					child2.SetGene(i, parent1.GetGene(j));
					break;
				}
				if (j == secondPoint - 1)
				{
					break;
				}
			}
			if (i == firstPoint - 1)
			{
				break;
			}
		}

		std::vector<Chromosome> children{ child1, child2 };
		for (auto& i : children)
		{
			i.CheckCorrectness();
		}
		return children;
	}

	static std::vector<Chromosome> CycleCrossover(const Chromosome& parent1, const Chromosome& parent2)
	{
		const std::size_t chromosomeLength = std::min(parent1.Length(), parent2.Length());
		Chromosome child1 = parent1;
		Chromosome child2 = parent2;

		std::vector<std::size_t> indices;
		std::size_t i = 0;
		indices.push_back(i);
		while (parent2.GetGene(i) != parent1.GetGene(0))
		{
			const std::vector<std::size_t> tmp = parent1.GetChromosome();
			i = std::distance(tmp.begin(), std::find(tmp.begin(), tmp.end(), parent2.GetGene(i)));
			indices.push_back(i);
		}
		for (std::size_t j = 0; j < chromosomeLength; ++j)
		{
			if (std::find(indices.begin(), indices.end(), j) != indices.end())
			{
				child1.SetGene(j, parent1.GetGene(j));
			}
			else
			{
				child1.SetGene(j, parent2.GetGene(j));
			}
		}

		indices.clear();
		i = 0;
		indices.push_back(i);

		while (parent1.GetGene(i) != parent2.GetGene(0))
		{
			const std::vector<std::size_t> tmp = parent2.GetChromosome();
			i = std::distance(tmp.begin(), std::find(tmp.begin(), tmp.end(), parent1.GetGene(i)));
			indices.push_back(i);
		}

		for (std::size_t j = 0; j < chromosomeLength; ++j)
		{
			if (std::find(indices.begin(), indices.end(), j) != indices.end())
			{
				child2.SetGene(j, parent2.GetGene(j));
			}
			else
			{
				child2.SetGene(j, parent1.GetGene(j));
			}
		}

		std::vector<Chromosome> children{ child1, child2 };
		for (auto& i : children)
		{
			i.CheckCorrectness();
		}
		return children;
	}
};
