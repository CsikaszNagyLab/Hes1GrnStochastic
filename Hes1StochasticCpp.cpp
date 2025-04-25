// Hes1StochasticCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <thread>
#include "CellModel.h"

#define NUM_CELLS 200

CellParameter initParameters() {
	CellParameter parameter;

	parameter.Volume = 1.3e-12;
	parameter.numTotalHes1Promoters = 15;
	parameter.numTotalProteasePromoters = 12;

	parameter.Hes1ProteinDecay = 0.108;
	parameter.hes1mRNADecay = 0.028;
	parameter.proteaseDecay = 0.1;
	parameter.dimerDecay = 0.036;

	parameter.hes1Transcription = 3;
	parameter.hes1Translation = 1;
	parameter.hes1TranscriptionReduction = 100;

	parameter.proteaseTranscription = 1000;
	parameter.proteaseTranscriptionReduction = 138;

	parameter.associationRate = 0.03;

	parameter.hes1PromoterBindRate = 0.5;
	parameter.hes1PromoterLiberationRate = 0.096;
	parameter.proteasePromoterBindRate = 150;
	parameter.proteasePromoterLiberationRate = 15;

	parameter.dimerizationRate = 1;

	parameter.bmalEffect = 0.7;

	return parameter;
}

struct simulationParameters {
	CellParameter parameter;
	CellState initialState;
	double simulationTime;
	double timeResolution;
};

struct simulationResult {
	std::vector<float> timeArray;
	std::vector<CellState> stateArray;

	void clear() {
		timeArray.clear();
		stateArray.clear();
	}
};

simulationResult RunSimulation(simulationParameters simParam, std::atomic<float>* progress) {
	float nextCheckPoint = simParam.timeResolution;
	CellModel model(simParam.parameter, simParam.initialState);

	std::vector<float> timeArray;
	timeArray.push_back(0);
	std::vector<CellState> stateArray;
	stateArray.push_back(simParam.initialState);

	do
	{
		model.Step();
		if (model.GetTimeElapsed() >= nextCheckPoint)
		{
			nextCheckPoint += simParam.timeResolution;
			//std::cout << "Time: " << std::fixed << std::setprecision(1) << model.GetTimeElapsed() << ", " << model.GetState().ToString() << std::endl;
			timeArray.push_back(model.GetTimeElapsed());
			stateArray.push_back(model.GetState());
			*progress = model.GetTimeElapsed() / simParam.simulationTime;
		}
	} while (model.GetTimeElapsed() <= simParam.simulationTime);
	//model.PrintReactionFireCounts();
	return { timeArray, stateArray };
}

void RunSimulationThread(simulationParameters simParam, simulationResult* result, std::atomic<float>* progress, std::atomic<bool>* done) {
	*result = RunSimulation(simParam, progress);
	*done = true;
}

int main()
{
	CellParameter parameter = initParameters();

	CellState initialState[NUM_CELLS];

	std::ifstream inFile("initialStates.csv");
	std::vector<CellState> cycleStates;
	std::string line;
	std::getline(inFile, line); // Skip header

	while (std::getline(inFile, line))
	{
		std::stringstream ss(line);
		std::string token;
		int stateValues[6];

		std::getline(ss, token, ',');
		float time = std::stof(token);

		for (int i = 0; i < 6; i++) {
			std::getline(ss, token, ',');
			stateValues[i] = std::stoi(token);
		}

		cycleStates.push_back({ 0, stateValues[0], stateValues[1], stateValues[2], stateValues[3], stateValues[4], stateValues[5] });
		ss.clear();
	}

	double simulationTime = 60*24*3;
	double timeResolution = .1;

	simulationResult cellResults[NUM_CELLS];
	std::vector<std::thread> threads;

	std::atomic<float> progresses[NUM_CELLS];
	std::atomic<bool> done[NUM_CELLS];

	const int SEED = 1;
	for (size_t i = 0; i < SEED; i++)
	{
		GetRandomNumber();
	}
	int percentages[3] = { 30 };
	//int percentages[3] = { 0, 70, 100 };
	for (int exp = 0; exp < sizeof(percentages) / sizeof(int); exp++)
	{
		parameter.bmalEffect = percentages[exp] / 100.0;

		for (size_t i = 0; i < NUM_CELLS; i++)
		{
			progresses[i] = 0;

			int stateId = round(GetRandomNumber() * (cycleStates.size() - 1));
			CellState initialState = cycleStates[stateId];
			simulationParameters simParam = { parameter, initialState, simulationTime, timeResolution };
			threads.push_back(std::thread(RunSimulationThread, simParam, &cellResults[i], &progresses[i], &done[i]));
		}

		bool allDone = false;
		do {
			for (int i = 0; i < NUM_CELLS; i++) {
				std::cout << "Cell " << std::setw(3) << i + 1 << ": " << std::fixed << std::setprecision(1) << std::setw(5) << progresses[i] * 100 << "%\t";
				if (i % 5 == 4) {
					std::cout << std::endl;
				}
			}
			std::cout << std::endl << "--------------------------------------------------------------------" << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(500));
			allDone = true;
			for (size_t i = 0; i < NUM_CELLS; i++)
			{
				allDone &= done[i];
			}
		} while (!allDone);

		for (auto& thread : threads) {
			thread.join();
		}

		for (int i = 0; i < NUM_CELLS; i++) {
			auto timeArray = cellResults[i].timeArray;
			auto stateArray = cellResults[i].stateArray;

			// Write to CSV file
			std::ofstream outFile("F_bmal" + std::to_string(percentages[exp]) + "_200/simulation_results_cell_" + std::to_string(i + 1) + ".csv");
			outFile << "Time,numHes1Protein,numHes1mRNA,numProtease,numHes1FreePromoter,numProteaseFreePromoter,numHes1Dimer\n";
			for (size_t i = 0; i < timeArray.size(); ++i) {
				outFile << timeArray[i] << ","
					<< stateArray[i].numHes1Protein << ","
					<< stateArray[i].numHes1mRNA << ","
					<< stateArray[i].numProtease << ","
					<< stateArray[i].numHes1FreePromoter << ","
					<< stateArray[i].numProteaseFreePromoter << ","
					<< stateArray[i].numHes1Dimer << "\n";
			}
			outFile.close();
			threads[i].~thread();
			cellResults[i].clear();
		}
		threads.clear();
		std::fill(done, done + NUM_CELLS, false);
		std::fill(progresses, progresses + NUM_CELLS, 0);
	}

	return 0;
}