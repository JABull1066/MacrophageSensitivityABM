/*

Copyright (c) 2005-2020, University of Oxford.
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

#include "MacrophagePhenotypeSwitchingCellCycle.hpp"
#include "RandomNumberGenerator.hpp"

MacrophagePhenotypeSwitchingCellCycle::MacrophagePhenotypeSwitchingCellCycle()
    : AbstractCellCycleModel(),
    mPhagocytosisCooldownPeriod(10.0),
    mTimeUntilPhagocytosisIsPossible(mPhagocytosisCooldownPeriod)
{
}

bool MacrophagePhenotypeSwitchingCellCycle::ReadyToDivide()
{
    // Update macrophage phenotype each timestep
    UpdatePhenotype();
    // If still in phagocytosis recharge period, reduce that timer
    if (mTimeUntilPhagocytosisIsPossible > 0)
    {
        mTimeUntilPhagocytosisIsPossible = mTimeUntilPhagocytosisIsPossible - SimulationTime::Instance()->GetTimeStep();
    }

    return false;
}

void MacrophagePhenotypeSwitchingCellCycle::UpdatePhenotype() {
    // Get cell's csf1 concentration
    double csf1_concentration = mpCell->GetCellData()->GetItem("csf1");
    // Get cell's phenotype
    double phenotype = mpCell->GetCellData()->GetItem("phenotype");

    if (csf1_concentration > 1)
    {
        phenotype = phenotype+0.0001;
    }
    else
    {
        phenotype = phenotype-0.0001;
    }
    // Add some drift
    double driftMean = 0;
    double driftStdDev = 0.001;
    phenotype = phenotype + RandomNumberGenerator::Instance()->NormalRandomDeviate(driftMean,driftStdDev);

    // Cap phenotype range
    phenotype = std::min(phenotype,1.0);
    phenotype = std::max(phenotype,0.0);
    mpCell->GetCellData()->SetItem("phenotype",phenotype);
}

// LCOV_EXCL_START
AbstractCellCycleModel* MacrophagePhenotypeSwitchingCellCycle::CreateCellCycleModel()
{
    NEVER_REACHED;
    return nullptr;
}
// LCOV_EXCL_STOP

double MacrophagePhenotypeSwitchingCellCycle::GetAverageTransitCellCycleTime()
{
    return DBL_MAX;
}

double MacrophagePhenotypeSwitchingCellCycle::GetAverageStemCellCycleTime()
{
    return DBL_MAX;
}


double MacrophagePhenotypeSwitchingCellCycle::GetPhagocytosisCooldownPeriod()
{
    return mPhagocytosisCooldownPeriod;
}

void MacrophagePhenotypeSwitchingCellCycle::SetPhagocytosisCooldownPeriod(double newPhagocytosisCooldownPeriod)
{
    mPhagocytosisCooldownPeriod = newPhagocytosisCooldownPeriod;
}

double MacrophagePhenotypeSwitchingCellCycle::GetTimeUntilPhagocytosisIsPossible()
{
    return mTimeUntilPhagocytosisIsPossible;
}

void MacrophagePhenotypeSwitchingCellCycle::SetTimeUntilPhagocytosisIsPossible(double newTimeUntilPhagocytosisIsPossible)
{
    mTimeUntilPhagocytosisIsPossible = newTimeUntilPhagocytosisIsPossible;
}


void MacrophagePhenotypeSwitchingCellCycle::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class

    *rParamsFile << "\t\t\t<PhagocytosisCooldownPeriod>" << mPhagocytosisCooldownPeriod << "</PhagocytosisCooldownPeriod>\n";

    // Call method on parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MacrophagePhenotypeSwitchingCellCycle)
