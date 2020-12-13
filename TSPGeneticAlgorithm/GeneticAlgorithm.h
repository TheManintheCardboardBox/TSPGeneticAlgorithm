#pragma once

#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "Breeding.h"
#include "Crossover.h"
#include "GenerationStrategy.h"
#include "Mutation.h"
#include "Population.h"
#include "Selection.h"
#include "Utils.h"

class GeneticAlgorithm
{
public:
	using DistanceMatrix = std::vector<std::vector<double>>;
	GeneticAlgorithm(const DistanceMatrix& distanceMatrix) :
		m_distanceMatrix(distanceMatrix),
		m_chromosomeLength(m_distanceMatrix.size())
	{
		SetInitialPopulationHeuristic();
		SetPopulationSize();
		SetBreedingSystem();
		SetInbreedingValue();
		SetAutobreedingValue();
		SetCrossoverOperator();
		SetCrossoverRate();
		SetMutationOperator();
		SetMutationRate();
		SetGenerationStrategy();
		SetOverlapRate();
		SetElitismUsage();
		SetSelectionOperator();
		SetTournamentSize();
		SetMaxGeneration();
		SetMaxPlateau();
	}

	void Run()
	{
		std::vector<std::string> inputString;
		std::string string;

		do
		{
			std::cout << "Available commands: run - run genetic algorithm" << std::endl;
			std::cout << "                    set - set a parameter" << std::endl;
			std::cout << "                    show - show genetic algorithm parameters" << std::endl;
			std::cout << "                    quit - quit the program" << std::endl;
			std::cout << "> ";

			std::getline(std::cin, string);
			std::istringstream iss(string);
			inputString.clear();
			while (iss >> string)
			{
				inputString.push_back(string);
			}

			switch (inputString.size())
			{
			case 1:
				{
					if (inputString.at(0) == "run")
					{
						std::cout << "Running..." << std::endl;
						RunAlgorithm();
					}
					else if (inputString.at(0) == "show")
					{
						PrintParameters();
					}
					else if (inputString.at(0) == "quit")
					{
						std::cout << "Quitting..." << std::endl;
					}
					else
					{
						std::cout << "Incorrect input" << std::endl << std::endl;
					}
				}
				break;
			case 4:
				{
					if (inputString.at(0) == "set")
					{
						if (inputString.at(1) == "initial" && inputString.at(2) == "population")
						{
							if (inputString.at(3) == "random")
							{
								SetInitialPopulationHeuristic(InitialPopulationHeuristic::random);
								std::cout << "Initial population heuristic is set to random generation" << std::endl << std::endl;
							}
							else if (inputString.at(3) == "NN")
							{
								SetInitialPopulationHeuristic(InitialPopulationHeuristic::nearestNeighbor);
								std::cout << "Initial population heuristic is set to nearest neighbor algorithm" << std::endl << std::endl;
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "population" && inputString.at(2) == "size")
						{
							if (IsNumber(inputString.at(3)))
							{
								const std::size_t generationSize = std::stoul(inputString.at(3));
								SetPopulationSize(generationSize);
								std::cout << "Population size is set to " << generationSize << std::endl << std::endl;
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "breeding" && inputString.at(2) == "system")
						{
							if (inputString.at(3) == "panmixia")
							{
								SetBreedingSystem(BreedingSystem::panmixia);
								std::cout << "Breeding system is set to " << inputString.at(3) << std::endl;
							}
							else if (inputString.at(3) == "inbreeding")
							{
								SetBreedingSystem(BreedingSystem::inbreeding);
								std::cout << "Breeding system is set to inbreeding" << std::endl;
								std::cout << "Set inbreeding value" << std::endl;
								std::cout << "> ";
								std::size_t inbreedingValue = 0;
								std::cin >> inbreedingValue;
								if (std::cin && inbreedingValue <= m_chromosomeLength)
								{
									SetInbreedingValue(inbreedingValue);
									std::cout << "Inbreeding value is set to " << inbreedingValue << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else if (inputString.at(3) == "autobreeding")
							{
								SetBreedingSystem(BreedingSystem::autobreeding);
								std::cout << "Breeding system is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set autobreeding value" << std::endl;
								std::cout << "> ";
								std::size_t autobreedingValue = 0;
								std::cin >> autobreedingValue;
								if (std::cin && autobreedingValue <= m_chromosomeLength)
								{
									SetAutobreedingValue(autobreedingValue);
									std::cout << "Autobreeding value is set to " << autobreedingValue << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "crossover" && inputString.at(2) == "operator")
						{
							if (inputString.at(3) == "OX")
							{
								SetCrossoverOperator(CrossoverOperator::OX);
								std::cout << "Crossover operator is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set crossover rate" << std::endl;
								std::cout << "> ";
								double crossoverRate = 0.0;
								std::cin >> crossoverRate;
								if (std::cin && crossoverRate >= 0.0 && crossoverRate <= 1.0)
								{
									SetCrossoverRate(crossoverRate);
									std::cout << "Crossover rate is set to "
										<< std::setprecision(5) << crossoverRate << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else if (inputString.at(3) == "CX")
							{
								SetCrossoverOperator(CrossoverOperator::CX);
								std::cout << "Crossover operator is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set crossover rate" << std::endl;
								std::cout << "> ";
								double crossoverRate = 0.0;
								std::cin >> crossoverRate;
								if (std::cin && crossoverRate >= 0.0 && crossoverRate <= 1.0)
								{
									SetCrossoverRate(crossoverRate);
									std::cout << "Crossover rate is set to "
										<< std::setprecision(5) << crossoverRate << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "mutation" && inputString.at(2) == "operator")
						{
							if (inputString.at(3) == "saltation")
							{
								SetMutationOperator(MutationOperator::saltation);
								std::cout << "Mutation operator is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set mutation rate" << std::endl;
								std::cout << "> ";
								double mutationRate = 0.0;
								std::cin >> mutationRate;
								if (std::cin && mutationRate >= 0.0 && mutationRate <= 1.0)
								{
									SetMutationRate(mutationRate);
									std::cout << "Mutation rate is set to "
										<< std::setprecision(5) << mutationRate << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else if (inputString.at(3) == "inversion")
							{
								SetMutationOperator(MutationOperator::inversion);
								std::cout << "Mutation operator is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set mutation rate" << std::endl;
								std::cout << "> ";
								double mutationRate = 0.0;
								std::cin >> mutationRate;
								if (std::cin && mutationRate >= 0.0 && mutationRate <= 1.0)
								{
									SetMutationRate(mutationRate);
									std::cout << "Mutation rate is set to "
										<< std::setprecision(5) << mutationRate << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "generation" && inputString.at(2) == "strategy")
						{
							if (inputString.at(3) == "non-overlapping")
							{
								SetGenerationStrategy(GenerationStrategy::nonOverlapping);
								std::cout << "Generation strategy is set to " << inputString.at(3) << std::endl;
							}
							else if (inputString.at(3) == "overlapping")
							{
								SetGenerationStrategy(GenerationStrategy::overlapping);
								std::cout << "Generation strategy is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set overlap rate" << std::endl;
								std::cout << "> ";
								double overlapRate = 0.0;
								std::cin >> overlapRate;
								if (std::cin && overlapRate >= 0.0 && overlapRate <= 1.0)
								{
									SetOverlapRate(overlapRate);
									std::cout << "Overlap rate is set to "
										<< std::setprecision(5) << overlapRate << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "elitism" && inputString.at(2) == "usage")
						{
							if (inputString.at(3) == "yes")
							{
								SetElitismUsage(true);
								std::cout << "Elitism usage is set to yes" << std::endl << std::endl;
							}
							else if (inputString.at(3) == "no")
							{
								SetElitismUsage(false);
								std::cout << "Elitism usage is set to no" << std::endl << std::endl;
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "generation" && inputString.at(2) == "strategy")
						{
							if (inputString.at(3) == "non-overlapping")
							{
								SetGenerationStrategy(GenerationStrategy::nonOverlapping);
								std::cout << "Generation strategy is set to " << inputString.at(3) << std::endl << std::endl;
							}
							else if (inputString.at(3) == "overlapping")
							{
								SetGenerationStrategy(GenerationStrategy::overlapping);
								std::cout << "Generation strategy is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set overlap rate" << std::endl;
								std::cout << "> ";
								double overlapRate = 0.0;
								std::cin >> overlapRate;
								if (std::cin && overlapRate >= 0.0 && overlapRate <= 1.0)
								{
									SetOverlapRate(overlapRate);
									std::cout << "Overlap rate is set to "
										<< std::setprecision(5) << overlapRate << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "selection" && inputString.at(2) == "operator")
						{
							if (inputString.at(3) == "proportionate")
							{
								SetSelectionOperator(SelectionOperator::fitnessProportionate);
								std::cout << "Selection operator is set to fitness proportionate selection" << std::endl << std::endl;
							}
							else if (inputString.at(3) == "tournament")
							{
								SetSelectionOperator(SelectionOperator::tournament);
								std::cout << "Selection operator is set to " << inputString.at(3) << std::endl << std::endl;
								std::cout << "Set tournament size" << std::endl;
								std::cout << "> ";
								std::size_t tournamentSize = 0;
								std::cin >> tournamentSize;
								if (std::cin && tournamentSize <= m_populationSize)
								{
									SetTournamentSize(tournamentSize);
									std::cout << "Tournament size is set to " << tournamentSize << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
								std::cin.ignore();
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "max" && inputString.at(2) == "generation")
						{
							if (IsNumber(inputString.at(3)))
							{
								const std::size_t maxGeneration = std::stoul(inputString.at(3));
								SetMaxGeneration(maxGeneration);
								std::cout << "Max generation is set to " << maxGeneration << std::endl << std::endl;
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else if (inputString.at(1) == "max" && inputString.at(2) == "plateau")
						{
							if (IsNumber(inputString.at(3)))
							{
								const std::size_t maxPlateau = std::stoul(inputString.at(3));
								if (std::cin && maxPlateau <= m_maxGeneration)
								{
									SetMaxPlateau(maxPlateau);
									std::cout << "Max plateau is set to " << maxPlateau << std::endl << std::endl;
								}
								else
								{
									std::cout << "Incorrect input" << std::endl << std::endl;
								}
							}
							else
							{
								std::cout << "Incorrect input" << std::endl << std::endl;
							}
						}
						else
						{
							std::cout << "Incorrect input" << std::endl << std::endl;
						}
					}
					else
					{
						std::cout << "Incorrect input" << std::endl << std::endl;
					}
				}
				break;
			default:
				std::cout << "Incorrect input" << std::endl << std::endl;
				break;
			}
		} while (inputString.at(0) != "quit");
	}

private:
	void SetInitialPopulationHeuristic(
		const InitialPopulationHeuristic& initialPopulationHeuristic = InitialPopulationHeuristic::random)
	{
		m_initialPopulationHeuristic = initialPopulationHeuristic;
	}

	void SetPopulationSize(const std::size_t& populationSize = 20)
	{
		m_populationSize = populationSize;
	}

	void SetBreedingSystem(const BreedingSystem& breedingSystem = BreedingSystem::panmixia)
	{
		m_breedingSystem = breedingSystem;
	}

	void SetInbreedingValue(const std::size_t& inbreedingValue = 1)
	{
		m_inbreedingValue = inbreedingValue;
	}

	void SetAutobreedingValue(const std::size_t& autobreedingValue = 1)
	{
		m_autobreedingValue = autobreedingValue;
	}

	void SetCrossoverOperator(const CrossoverOperator& crossoverOperator = CrossoverOperator::OX)
	{
		m_crossoverOperator = crossoverOperator;
	}

	void SetCrossoverRate(const double& crossoverRate = 0.95)
	{
		m_crossoverRate = crossoverRate;
	}

	void SetMutationOperator(const MutationOperator& mutationOperator = MutationOperator::saltation)
	{
		m_mutationOperator = mutationOperator;
	}

	void SetMutationRate(const double& mutationRate = 0.05)
	{
		m_mutationRate = mutationRate;
	}

	void SetGenerationStrategy(const GenerationStrategy& generationStrategy = GenerationStrategy::nonOverlapping)
	{
		m_generationStrategy = generationStrategy;
	}

	void SetOverlapRate(const double& overlapRate = 0.5)
	{
		m_overlapRate = overlapRate;
	}

	void SetElitismUsage(const bool& isElitismUsed = true)
	{
		m_isElitismUsed = isElitismUsed;
	}

	void SetSelectionOperator(const SelectionOperator& selectionOperator = SelectionOperator::fitnessProportionate)
	{
		m_selectionOperator = selectionOperator;
	}

	void SetTournamentSize(const std::size_t& tournamentSize = 10)
	{
		m_tournamentSize = tournamentSize;
	}

	void SetMaxGeneration(const std::size_t& maxGeneration = 50)
	{
		m_maxGeneration = maxGeneration;
	}

	void SetMaxPlateau(const std::size_t& maxPlateau = 10)
	{
		m_maxPlateau = maxPlateau;
	}

	void RunAlgorithm()
	{
		Population initialPopulation(m_populationSize, m_distanceMatrix, m_initialPopulationHeuristic);
		const Chromosome initialPopulationFittest = initialPopulation.GetFittest();
		Chromosome currentBestSolution = initialPopulationFittest;
		PrintResults(0, initialPopulation);
		std::cout << "New current best solution: [ " << currentBestSolution << " ] - "
			<< std::setprecision(5) << currentBestSolution.GetDistance() << std::endl;

		Population nextGeneration;
		std::size_t plateauCount = 0;
		for (std::size_t i = 0; i < m_maxGeneration; ++i)
		{
			nextGeneration = NextGeneration(initialPopulation);
			PrintResults(i + 1, nextGeneration);
			const Chromosome nextGenerationFittest = nextGeneration.GetFittest();
			if (currentBestSolution.GetFitness() >= nextGenerationFittest.GetFitness())
			{
				plateauCount++;
				std::cout << "Current best solution: [ " << currentBestSolution << " ] - "
					<< std::setprecision(5) << currentBestSolution.GetDistance() << std::endl;
			}
			else
			{
				currentBestSolution = nextGenerationFittest;
				std::cout << "New current best solution: [ " << currentBestSolution << " ] - "
					<< std::setprecision(5) << currentBestSolution.GetDistance() << std::endl;
			}
			if (plateauCount == m_maxPlateau)
			{
				break;
			}
		}

		std::cout << "---------------------------------------------------------------------------" << std::endl;
		if (plateauCount == m_maxPlateau)
		{
			std::cout << "Algorithm reached a plateau" << std::endl;
		}
		else
		{
			std::cout << "Algorithm reached a maximum generation" << std::endl;
		}
		std::cout << "Solution: [ " << currentBestSolution << " ] - "
			<< std::setprecision(5) << currentBestSolution.GetDistance() << std::endl << std::endl;
	}

	Population NextGeneration(const Population& population)
	{
		Population reproductiveSet;
		for (std::size_t i = 0; i < m_populationSize; ++i)
		{
			std::vector<Chromosome> parents;
			switch (m_breedingSystem)
			{
			case BreedingSystem::panmixia:
				parents = Breeding::Panmixia(population);
				break;
			case BreedingSystem::inbreeding:
				parents = Breeding::Inbreeding(population, m_inbreedingValue);
				break;
			case BreedingSystem::autobreeding:
				parents = Breeding::Autobreeding(population, m_autobreedingValue);
				break;
			default:
				break;
			}
			const Chromosome parent1 = parents.at(0);
			const Chromosome parent2 = parents.at(1);

			std::vector<Chromosome> children;
			if (GenerateRandomReal(0.0, 1.0) <= m_crossoverRate)
			{
				switch (m_crossoverOperator)
				{
				case CrossoverOperator::OX:
					children = Crossover::OrderCrossover(parent1, parent2);
					break;
				case CrossoverOperator::CX:
					children = Crossover::CycleCrossover(parent1, parent2);
					break;
				default:
					break;
				}
			}
			else
			{
				children = parents;
			}
			for (auto& i : children)
			{
				if (GenerateRandomReal(0.0, 1.0) <= m_mutationRate)
				{
					switch (m_mutationOperator)
					{
					case MutationOperator::saltation:
						Mutation::SaltationMutation(i);
						break;
					case MutationOperator::inversion:
						Mutation::InversionMutation(i);
						break;
					default:
						break;
					}
				}
			}
			for (auto& i : children)
			{
				reproductiveSet.Add(i);
			}
		}

		Population nextGeneration;
		std::size_t elitismValue = 0;
		if (m_isElitismUsed)
		{
			nextGeneration.Add(population.GetFittest());
			elitismValue = 1;
		}
		switch (m_generationStrategy)
		{
		case GenerationStrategy::nonOverlapping:
			{
				for (std::size_t i = 0; i < m_populationSize - elitismValue; ++i)
				{
					switch (m_selectionOperator)
					{
					case SelectionOperator::fitnessProportionate:
						{
							const Chromosome candidate = Selection::FitnessProportionateSelection(reproductiveSet);
							nextGeneration.Add(candidate);
						}
						break;
					case SelectionOperator::tournament:
						{
							const Chromosome candidate
								= Selection::TournamentSelection(reproductiveSet, m_tournamentSize);
							nextGeneration.Add(candidate);
						}
						break;
					default:
						break;
					}
				}
			}
			break;
		case GenerationStrategy::overlapping:
			{
				const std::size_t reproductiveSize = static_cast<std::size_t>(m_overlapRate * m_populationSize);
				for (std::size_t i = 0; i < m_populationSize - reproductiveSize - elitismValue; ++i)
				{
					nextGeneration.Add(population.GetChromosome(GenerateRandomInt(0, m_populationSize - 1)));
				}

				for (std::size_t i = 0; i < reproductiveSize; ++i)
				{
					switch (m_selectionOperator)
					{
					case SelectionOperator::fitnessProportionate:
						{
							const Chromosome candidate = Selection::FitnessProportionateSelection(reproductiveSet);
							nextGeneration.Add(candidate);
						}
						break;
					case SelectionOperator::tournament:
						{
							const Chromosome candidate
								= Selection::TournamentSelection(reproductiveSet, m_tournamentSize);
							nextGeneration.Add(candidate);
						}
						break;
					default:
						break;
					}
				}
			}
			break;
		default:
			break;
		}
		return nextGeneration;
	}

	void PrintParameters()
	{
		std::cout << "Current algorithm parameters" << std::endl;

		std::cout << "Initial population heuristic: ";
		switch (m_initialPopulationHeuristic)
		{
		case InitialPopulationHeuristic::random:
			std::cout << "random generation" << std::endl;
			break;
		case InitialPopulationHeuristic::nearestNeighbor:
			std::cout << "nearest neighbor algorithm" << std::endl;
			break;
		default:
			break;
		}

		std::cout << "Population size: " << m_populationSize << std::endl;

		std::cout << "Breeding system: ";
		switch (m_breedingSystem)
		{
		case BreedingSystem::panmixia:
			std::cout << "panmixia" << std::endl;
			break;
		case BreedingSystem::inbreeding:
			std::cout << "inbreeding with value " << m_inbreedingValue << std::endl;
			break;
		case BreedingSystem::autobreeding:
			std::cout << "autobreeding with value " << m_autobreedingValue << std::endl;
			break;
		default:
			break;
		}

		std::cout << "Crossover operator: ";
		switch (m_crossoverOperator)
		{
		case CrossoverOperator::OX:
			std::cout << "OX with rate " << m_crossoverRate << std::endl;
			break;
		case CrossoverOperator::CX:
			std::cout << "CX with rate " << m_crossoverRate << std::endl;
			break;
		default:
			break;
		}

		std::cout << "Mutation operator: ";
		switch (m_mutationOperator)
		{
		case MutationOperator::saltation:
			std::cout << "saltation with rate " << m_mutationRate << std::endl;
			break;
		case MutationOperator::inversion:
			std::cout << "inversion with rate " << m_mutationRate << std::endl;
			break;
		default:
			break;
		}

		std::cout << "Generation strategy: ";
		switch (m_generationStrategy)
		{
		case GenerationStrategy::nonOverlapping:
			std::cout << "non-overlapping" << std::endl;
			break;
		case GenerationStrategy::overlapping:
			std::cout << "overlapping with rate " << m_overlapRate << std::endl;
			break;
		default:
			break;
		}
		std::cout << "Elitism usage: " << (m_isElitismUsed ? "yes" : "no") << std::endl;

		std::cout << "Selection operator: ";
		switch (m_selectionOperator)
		{
		case SelectionOperator::fitnessProportionate:
			std::cout << "fitness proportionate selection" << std::endl;
			break;
		case SelectionOperator::tournament:
			std::cout << "tournament with size " << m_tournamentSize << std::endl;
			break;
		default:
			break;
		}

		std::cout << "Max generation: " << m_maxGeneration << std::endl;
		std::cout << "Max plateau generation: " << m_maxPlateau << std::endl << std::endl;
	}

	void PrintResults(const std::size_t& numGeneration, Population& generation)
	{
		std::cout << "---------------------------------------------------------------------------" << std::endl;
		std::cout << "Generation " << numGeneration << ": ";

		std::string delim = "";
		const std::size_t numDigit = (numGeneration >= 10) ? (numGeneration >= 100) ? 3 : 2 : 1;
		const std::string indent = std::string((11 + numDigit + 2), ' ');

		generation.Shuffle();
		for (auto& i : generation.GetPopulation())
		{
			std::cout << delim << "[ " << i << " ] - " << std::setprecision(5) << i.GetDistance() << std::endl;
			delim = indent;
		}

		const Chromosome currentBestSolution = generation.GetFittest();
		std::cout << "Fittest of the generation: [ " << currentBestSolution << " ] - "
			<< std::setprecision(5) << currentBestSolution.GetDistance() << std::endl;
	}

	InitialPopulationHeuristic m_initialPopulationHeuristic;
	std::size_t m_populationSize;
	BreedingSystem m_breedingSystem;
	std::size_t m_inbreedingValue;
	std::size_t m_autobreedingValue;
	CrossoverOperator m_crossoverOperator;
	double m_crossoverRate;
	MutationOperator m_mutationOperator;
	double m_mutationRate;
	GenerationStrategy m_generationStrategy;
	double m_overlapRate;
	bool m_isElitismUsed;
	SelectionOperator m_selectionOperator;
	std::size_t m_tournamentSize;
	std::size_t m_maxGeneration;
	std::size_t m_maxPlateau;

	std::vector<std::vector<double>> m_distanceMatrix;
	std::size_t m_chromosomeLength;
};
