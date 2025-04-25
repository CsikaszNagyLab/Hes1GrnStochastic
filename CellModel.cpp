#include "CellModel.h"
#include <math.h>

CellModel::CellModel(CellParameter parameter, CellState initialState)
{
	this->parameter = parameter;
	this->state = initialState;
	this->state.timeElapsed = 0;
	CalculateAllPropensities();

	for (size_t i = 0; i < NUM_REACTIONS; i++)
	{
		this->reactionFireCounts[i] = 0;
	}
}

CellModel::~CellModel()
{

}

void CellModel::SetState(CellState state)
{
	this->state = state;
}

CellState CellModel::GetState()
{
	return this->state;
}

long double CellModel::GetTimeElapsed()
{
	return this->state.timeElapsed;
}

void CellModel::PrintReactionFireCounts()
{
	for (int i = 0; i < NUM_REACTIONS; i++)
	{
		printf("Recation %d: %d\n", i + 1, this->reactionFireCounts[i]);
	}
}

void CellModel::CalculateAllPropensities()
{
	RecalculatePropensities();
}

void CellModel::ExecuteSingleReaction(int reactionId)
{
	switch (reactionId)
	{
		case 0:
			state.numHes1mRNA--;
			break;
		case 1:
			state.numHes1Protein--;
			break;
		case 2:
			state.numHes1Dimer--;
			break;
		case 3:
			state.numProtease--;
			break;
		case 4:
			state.numHes1Protein--;
			state.numProtease--;
			break;
		case 5:
			state.numHes1Protein++;
			break;
		case 6:
			state.numHes1Protein -= 2;
			state.numHes1Dimer++;
			break;
		case 7:
			state.numHes1mRNA++;
			break;
		case 8:
			state.numHes1mRNA++;
			break;
		case 9:
			state.numProtease++;
			break;
		case 10:
			state.numProtease++;
			break;
		case 11:
			state.numHes1FreePromoter--;
			state.numHes1Dimer--;
			break;
		case 12:
			state.numHes1FreePromoter++;
			break;
		case 13:
			state.numProteaseFreePromoter--;
			state.numHes1Dimer--;
			break;
		case 14:
			state.numProteaseFreePromoter++;
			break;
	}
	RecalculatePropensities();
	if (state.anyNegative())
	{
		throw std::runtime_error("Negative state");
	}
}

void CellModel::ExecuteReaction()
{
	double sum = 0;
	for (int i = 0; i < 15; i++)
	{
		sum += propensities[i];
	}
	double r1 = GetRandomNumber();
	double r2 = GetRandomNumber();
	long double tau = static_cast<long double>(-log(r1)) / sum;

	state.timeElapsed += tau;

	double sum2 = 0;
	for (int i = 0; i < 15; i++)
	{
		sum2 += propensities[i];
		if (r2 * sum <= sum2)
		{
			ExecuteSingleReaction(i);
			this->reactionFireCounts[i]++;
			break;
		}
	}
}

void CellModel::Step() {
	ExecuteReaction();
}

void CellModel::RecalculatePropensities() {
	float K = parameter.Volume * 6e14;
	double bmalLevel = sin(2 * 3.14159265 * (state.timeElapsed + parameter.cycleOffset) / (60 * 24)) / 2 + 0.5;
	double bmalEffect = (1 + parameter.bmalEffect * bmalLevel);

	propensities[7] = parameter.hes1Transcription * state.numHes1FreePromoter * bmalEffect;
	propensities[8] = parameter.hes1Transcription / parameter.hes1TranscriptionReduction
		* (parameter.numTotalHes1Promoters - state.numHes1FreePromoter) * bmalEffect;

	propensities[0] = parameter.hes1mRNADecay * state.numHes1mRNA;
	propensities[1] = parameter.Hes1ProteinDecay * state.numHes1Protein;
		
	propensities[2] = parameter.dimerDecay * state.numHes1Dimer;
		
	propensities[3] = parameter.proteaseDecay * state.numProtease;
		
	propensities[4] = parameter.associationRate * state.numHes1Protein * state.numProtease / K;
		
	propensities[5] = parameter.hes1Translation * state.numHes1mRNA;
		
	propensities[6] = parameter.dimerizationRate * state.numHes1Protein * (state.numHes1Protein - 1) / K;
		
	propensities[9] = parameter.proteaseTranscription * state.numProteaseFreePromoter;
	propensities[10] = parameter.proteaseTranscription / parameter.proteaseTranscriptionReduction
		* (parameter.numTotalProteasePromoters - state.numProteaseFreePromoter);
		propensities[11] = parameter.hes1PromoterBindRate * state.numHes1Dimer * state.numHes1FreePromoter / K;
	propensities[12] = parameter.hes1PromoterLiberationRate * (parameter.numTotalHes1Promoters - state.numHes1FreePromoter);
	propensities[13] = parameter.proteasePromoterBindRate * state.numHes1Dimer * state.numProteaseFreePromoter / K;
		propensities[14] = parameter.proteasePromoterLiberationRate * (parameter.numTotalProteasePromoters - state.numProteaseFreePromoter);
	int i = 0;
}