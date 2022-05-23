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

#include "Arwert2018_MacrophageCellCycle.hpp"
#include "RandomNumberGenerator.hpp"

Arwert2018_MacrophageCellCycle::Arwert2018_MacrophageCellCycle()
    : AbstractCellCycleModel(),
    mPhagocytosisCooldownPeriod(10.0),
    mCriticalTGFbetaThreshold(2),
    mPhenotypeIncrementPerHour(0.1),
    mTimeUntilPhagocytosisIsPossible(mPhagocytosisCooldownPeriod)
{
}

bool Arwert2018_MacrophageCellCycle::ReadyToDivide()
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

void Arwert2018_MacrophageCellCycle::UpdatePhenotype() {
    // Get cell's csf1 concentration
    double tgf_concentration = mpCell->GetCellData()->GetItem("tgf");
    // Get cell's phenotype
    double phenotype = mpCell->GetCellData()->GetItem("phenotype");

    // Phenotype increases by a small amount each time step
    // todo Update this with something more realistic

    if (tgf_concentration > mCriticalTGFbetaThreshold)
    {
        double dt = SimulationTime::Instance()->GetTimeStep();
        phenotype = phenotype + mPhenotypeIncrementPerHour*dt;
    }

//    double egf;
//    if (phenotype > 0.5)
//    {
//        // If macrophage is M2, it produces EGF
//        egf = mpCell->GetCellData()->GetItem("egf");
//        egf = egf + 0.0001;
//    }
//    else{
//        // If macrophage is M1, EGF decays
//        egf = mpCell->GetCellData()->GetItem("egf");
//        egf = egf - 0.0001;
//    }

    // Cap phenotype range
    phenotype = std::min(phenotype,1.0);
    phenotype = std::max(phenotype,0.0);
    mpCell->GetCellData()->SetItem("phenotype",phenotype);

//    // Cap egf range
//    egf = std::min(egf,1.0);
//    egf = std::max(egf,0.0);
//    mpCell->GetCellData()->SetItem("egf",egf);

}

// LCOV_EXCL_START
AbstractCellCycleModel* Arwert2018_MacrophageCellCycle::CreateCellCycleModel()
{
    NEVER_REACHED;
    return nullptr;
}
// LCOV_EXCL_STOP

double Arwert2018_MacrophageCellCycle::GetAverageTransitCellCycleTime()
{
    return DBL_MAX;
}

double Arwert2018_MacrophageCellCycle::GetAverageStemCellCycleTime()
{
    return DBL_MAX;
}


double Arwert2018_MacrophageCellCycle::GetCriticalTGFbetaThreshold()
{
    return mCriticalTGFbetaThreshold;
}

void Arwert2018_MacrophageCellCycle::SetCriticalTGFbetaThreshold(double newCriticalTGFbetaThreshold)
{
    mCriticalTGFbetaThreshold = newCriticalTGFbetaThreshold;
}

double Arwert2018_MacrophageCellCycle::GetPhenotypeIncrementPerHour()
{
    return mPhenotypeIncrementPerHour;
}

void Arwert2018_MacrophageCellCycle::SetPhenotypeIncrementPerHour(double newPhenotypeIncrementPerHour)
{
    mPhenotypeIncrementPerHour = newPhenotypeIncrementPerHour;
}
double Arwert2018_MacrophageCellCycle::GetPhagocytosisCooldownPeriod()
{
    return mPhagocytosisCooldownPeriod;
}

void Arwert2018_MacrophageCellCycle::SetPhagocytosisCooldownPeriod(double newPhagocytosisCooldownPeriod)
{
    mPhagocytosisCooldownPeriod = newPhagocytosisCooldownPeriod;
}

double Arwert2018_MacrophageCellCycle::GetTimeUntilPhagocytosisIsPossible()
{
    return mTimeUntilPhagocytosisIsPossible;
}

void Arwert2018_MacrophageCellCycle::SetTimeUntilPhagocytosisIsPossible(double newTimeUntilPhagocytosisIsPossible)
{
    mTimeUntilPhagocytosisIsPossible = newTimeUntilPhagocytosisIsPossible;
}


void Arwert2018_MacrophageCellCycle::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class

    *rParamsFile << "\t\t\t<PhagocytosisCooldownPeriod>" << mPhagocytosisCooldownPeriod << "</PhagocytosisCooldownPeriod>\n";
    *rParamsFile << "\t\t\t<mPhenotypeIncrementPerHour>" << mPhenotypeIncrementPerHour << "</mPhenotypeIncrementPerHour>\n";
    *rParamsFile << "\t\t\t<mCriticalTGFbetaThreshold>" << mCriticalTGFbetaThreshold << "</mCriticalTGFbetaThreshold>\n";

    // Call method on parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Arwert2018_MacrophageCellCycle)
