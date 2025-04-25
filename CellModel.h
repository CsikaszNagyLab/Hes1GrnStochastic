#pragma once

#include <random>
#include <mutex>
#include <sstream>
#include <string>
#include <iomanip>

static std::mt19937 randomEngine;
static std::mutex randomMutex;

static double GetRandomNumber()
{
	std::lock_guard<std::mutex> lock(randomMutex);
	std::uniform_real_distribution<double> distribution(std::numeric_limits<double>::min(), 1);
	return distribution(randomEngine);
}

struct CellState
{
	long double timeElapsed;

	int numHes1Protein;
	int numHes1mRNA;
	int numProtease;
	int numHes1FreePromoter;
	int numProteaseFreePromoter;
	int numHes1Dimer;

	std::string ToString() const {
		std::ostringstream oss;
		oss << "State: "
			<< "Hes1 Prot.: " << std::setw(5) << numHes1Protein << ", "
			<< "hes1 mRNA: " << std::setw(5) << numHes1mRNA << ", "
			<< "Protease: " << std::setw(5) << numProtease << ", "
			<< "Free hes1 prom.: " << std::setw(3) << numHes1FreePromoter << ", "
			<< "Free prot. prom.: " << std::setw(3) << numProteaseFreePromoter << ", "
			<< "Hes1 Dimer: " << std::setw(5) << numHes1Dimer
			<< " }";
		return oss.str();
	}

	bool anyNegative() const {
		return numHes1Protein < 0 || numHes1mRNA < 0 || numProtease < 0 || numHes1FreePromoter < 0 || numProteaseFreePromoter < 0 || numHes1Dimer < 0;
	};
};

struct CellParameter {
	double Volume;

	int numTotalHes1Promoters;
	int numTotalProteasePromoters;

	double Hes1ProteinDecay;
	double hes1mRNADecay;
	double proteaseDecay;
	double dimerDecay;

	double hes1Transcription;
	double hes1Translation;
	double hes1TranscriptionReduction;

	double proteaseTranscription;
	double proteaseTranscriptionReduction;

	double associationRate;

	double proteasePromoterBindRate;
	double hes1PromoterBindRate;
	double proteasePromoterLiberationRate;
	double hes1PromoterLiberationRate;

	double dimerizationRate;
	double bmalEffect;

	double cycleOffset;
};

class CellModel
{
#define NUM_REACTIONS 15
private:
	CellState state;
	CellParameter parameter;
	double propensities[NUM_REACTIONS];

	void CalculateAllPropensities();
	void RecalculatePropensities();
	void ExecuteReaction();
	void ExecuteSingleReaction(int reactionId);

	int reactionFireCounts[NUM_REACTIONS];
public:
	CellModel(CellParameter parameter, CellState initialState);
	~CellModel();

	void SetState(CellState state);
	CellState GetState();
	long double GetTimeElapsed();

	void PrintReactionFireCounts();

	void Step();
};