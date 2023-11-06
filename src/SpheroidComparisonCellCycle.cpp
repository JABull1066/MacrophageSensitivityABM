/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "SpheroidComparisonCellCycle.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

#include "ApoptoticCellProperty.hpp"
#include "Debug.hpp"

SpheroidComparisonCellCycle::SpheroidComparisonCellCycle()
        : AbstractSimpleCellCycleModel(),
          mMinCellCycleDuration(12.0), // Hours
          mMaxCellCycleDuration(14.0),  // Hours
          mCellCycleClock(0),
          mCellCycleClockIncrement(1/120),
          mHypoxicConcentration(0.5),
          mNecroticConcentration(0.3)
{

}


SpheroidComparisonCellCycle::SpheroidComparisonCellCycle(const SpheroidComparisonCellCycle& rModel)
        : AbstractSimpleCellCycleModel(rModel),
          mMinCellCycleDuration(rModel.mMinCellCycleDuration),
          mMaxCellCycleDuration(rModel.mMaxCellCycleDuration),
          mCellCycleClock(rModel.mCellCycleClock),
          mCellCycleClockIncrement(rModel.mCellCycleClockIncrement),
          mHypoxicConcentration(rModel.mHypoxicConcentration),
          mNecroticConcentration(rModel.mNecroticConcentration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* SpheroidComparisonCellCycle::CreateCellCycleModel()
{
    return new SpheroidComparisonCellCycle(*this);
}

void SpheroidComparisonCellCycle::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCellCycleDuration = DBL_MAX;
        mCellCycleClockIncrement = 0;
    }
    else
    {
        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf(); // U[MinCCD,MaxCCD]

        // Ensure cell cycle clock is calibrated initially to age of cell (cells don't start sims at age 0)
        mCellCycleClock = GetAge() / mCellCycleDuration;
        // Work out increment clock should increase by each timestep, ensuring it hits 1 once mCellCycleDuration is reached
        //mCellCycleClockIncrement = SimulationTime::Instance()->GetTimeStep() / mCellCycleDuration;
        //TODO Assumes timestep is 1/120
        mCellCycleClockIncrement = 1 / (120*mCellCycleDuration);

    }

}

double SpheroidComparisonCellCycle::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void SpheroidComparisonCellCycle::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double SpheroidComparisonCellCycle::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void SpheroidComparisonCellCycle::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double SpheroidComparisonCellCycle::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double SpheroidComparisonCellCycle::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}


double SpheroidComparisonCellCycle::GetHypoxicConcentration() const
{
    return mHypoxicConcentration;
}

void SpheroidComparisonCellCycle::SetHypoxicConcentration(double hypoxicConcentration)
{
    assert(hypoxicConcentration<=1.0);
    assert(hypoxicConcentration>=0.0);
    mHypoxicConcentration = hypoxicConcentration;
}

double SpheroidComparisonCellCycle::GetNecroticConcentration() const
{
    return mNecroticConcentration;
}

void SpheroidComparisonCellCycle::SetNecroticConcentration(double necroticConcentration)
{
    assert(necroticConcentration <= 1.0);
    assert(necroticConcentration >= 0.0);
    mNecroticConcentration = necroticConcentration;
}

void SpheroidComparisonCellCycle::UpdateHypoxicCellCycles()
{
    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    // Get cell's oxygen concentration
    double oxygen_concentration = mpCell->GetCellData()->GetItem("oxygen");

    // Increment cell cycle clock by relevant proportion
    mCellCycleClock = mCellCycleClock + mCellCycleClockIncrement;

    // If cell is hypoxic, freeze its cell cycle. We can do this by either decreasing its age or increasing its cell cycle duration. Here, we increase the cell cycle duration
    if (oxygen_concentration < mHypoxicConcentration)
    {
        // Increase cell cycle duration by length of timestep
        mCellCycleDuration = mCellCycleDuration + SimulationTime::Instance()->GetTimeStep();

        // Clock does not increment
        mCellCycleClock = mCellCycleClock - mCellCycleClockIncrement;
    }
    // If cell is necrotic, kill it
    if (oxygen_concentration < mNecroticConcentration) {
        boost::shared_ptr <AbstractCellProperty> p_apoptotic_property =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        mpCell->AddCellProperty(p_apoptotic_property);
    }
}

// Overridden ReadyToDivide
bool SpheroidComparisonCellCycle::ReadyToDivide()
{
    assert(mpCell != nullptr);
    UpdateHypoxicCellCycles();
    if (!mReadyToDivide)
    {
        if (GetAge() >= mCellCycleDuration )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}


double SpheroidComparisonCellCycle::GetCellCycleProgression()
{
//    MARK;
//    PRINT_VARIABLE(mCellCycleClock);
//    PRINT_VARIABLE(mCellCycleClockIncrement);
    return mCellCycleClock;
}

//double SpheroidComparisonCellCycle::GetAge()
//{
//    MARK;
//    PRINT_VARIABLE(mCellCycleClock);
//    PRINT_VARIABLE(mCellCycleClockIncrement);
//    return mCellCycleClock;
//}

void SpheroidComparisonCellCycle::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<HypoxicConcentration>" << mHypoxicConcentration << "</HypoxicConcentration>\n";
    *rParamsFile << "\t\t\t<NecroticConcentration>" << mNecroticConcentration << "</NecroticConcentration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SpheroidComparisonCellCycle)
