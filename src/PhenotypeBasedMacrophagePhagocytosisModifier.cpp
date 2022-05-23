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

#include "PhenotypeBasedMacrophagePhagocytosisModifier.hpp"
#include "MacrophagePhenotypeSwitchingCellCycle.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellLabel.hpp"
#include "SimulationTime.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::PhenotypeBasedMacrophagePhagocytosisModifier()
    : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>(),
    mMacrophageLabel(7),
    mCellTypesToEat({0})
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::~PhenotypeBasedMacrophagePhagocytosisModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{

    auto nbcp = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
//    dynamic_cast <PottsBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation()));
    /*
     * Loop over the cell population to find macrophages.
     * At each time step, each macrophage checks a neighbouring cell at determines whether to begin phagocytosis
     */
    // Iterate over cell population
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Is it a macrophage?
        unsigned label = 0;
        CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection();

        if (collection.HasProperty<CellLabel>())
        {
            CellPropertyCollection collectionLabel = collection.GetProperties<CellLabel>();
//            collectionLabel.GetProperty()->GetColour();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collectionLabel.GetProperty());
            label = p_label->GetColour();


            if (mMacrophageLabel == label) {


                auto cyclemodel = static_cast<MacrophagePhenotypeSwitchingCellCycle*>((*cell_iter)->GetCellCycleModel());
                bool canEat = cyclemodel->GetTimeUntilPhagocytosisIsPossible() <= 0;
                // Do stuff

                if (canEat)
                {
                    // Determine whether eating will happen based on phenotype
                    // Use phenomenological hill function for prob per hour:
                    // half maximal at 0.65 - phenotype less than 0.5 means phagocytosis is likely
                    // Phenotype close to 1 means phagocytosis is unlikely
                    /*
                     * We assume a constant time step and that there are an integer number (n = 1/dt)
                     * of time steps per hour. We also assume that this method is called every time step
                     * and that the probabilities of dying at different times are independent.
                     *
                     * Let q=prob of eating per hour and p="probability of eating in a given time step".
                     *
                     * Probability of not eating in an hour:
                     * (1-q) = (1-p)^n = (1-p)^(1/dt).
                     *
                     * Rearranging for p:
                     * p = 1 - (1-q)^dt.
                     */
                    double x = (*cell_iter)->GetCellData()->GetItem("phenotype");
//                    double probOfEatingPerHour = 1 - (pow(x,10) / (pow(0.65,10) + pow(x,10)));
                    double maxProbOfEatingPerHour = 0.2;

                    double probOfEatingPerHour = maxProbOfEatingPerHour * (1 - ( pow(x,10) / ( pow(x,10) + pow(0.5,10) ) ));
                    double eating_prob_this_timestep = 1.0 - pow((1.0 - probOfEatingPerHour), SimulationTime::Instance()->GetTimeStep());

                    if (RandomNumberGenerator::Instance()->ranf() < eating_prob_this_timestep)
                    {
                        // Great, let's find some macrophage food
                        unsigned index = (&rCellPopulation)->GetLocationIndexUsingCell(*cell_iter);
                        // Neighbours for phagocytosis are within distance 1
                        std::set<unsigned> neighbours = nbcp->GetNodesWithinNeighbourhoodRadius(index, 1);
                        std::vector<unsigned> v(neighbours.begin(), neighbours.end());

                        std::random_shuffle ( v.begin(), v.end() );
                        // Now check each shuffled index in turn - when we find a cell we can kill, we eat it
                        for (unsigned i = 0; i < v.size(); i++)
                        {
                            double candidateIndex = v[i];
                            CellPtr p_targetCell = rCellPopulation.GetCellUsingLocationIndex(candidateIndex);
                            // Is cell a living tumour cell?
                            CellPropertyCollection collectionLabelTarget = p_targetCell->rGetCellPropertyCollection().GetProperties<CellLabel>();
                            boost::shared_ptr<CellLabel> p_labelTarget = boost::static_pointer_cast<CellLabel>(collectionLabelTarget.GetProperty());
                            unsigned targetlabel = p_labelTarget->GetColour();

                            bool isTargetLabelEdible = std::find(mCellTypesToEat.begin(), mCellTypesToEat.end(), targetlabel) != mCellTypesToEat.end();
                            if (!p_targetCell->HasApoptosisBegun() && isTargetLabelEdible)//targetlabel == 3)
                            {
                                p_targetCell->StartApoptosis();
                                // Reset phagocytosis timer
                                cyclemodel->SetTimeUntilPhagocytosisIsPossible(cyclemodel->GetPhagocytosisCooldownPeriod());
                                i = v.size();
                                continue;
                            }

                        }
                    }
                }
                else
                {
                    // Update phagocytosis cooldown period
                    cyclemodel->SetTimeUntilPhagocytosisIsPossible( cyclemodel->GetTimeUntilPhagocytosisIsPossible() - SimulationTime::Instance()->GetTimeStep());
                }

            }
        }

    }

    UpdateCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
}


//    // Iterate over cell population
//    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
//         cell_iter != rCellPopulation.End();
//         ++cell_iter)
//    {
//        // Get the volume of this cell
//        double cell_volume = rCellPopulation.GetVolumeOfCell(*cell_iter);
//
//        // Store the cell's volume in CellData
//        cell_iter->GetCellData()->SetItem("volume", cell_volume);
//    }

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::GetMacrophageLabel()
{
    return mMacrophageLabel;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::SetMacrophageLabel(unsigned newMacrophageLabel)
{
    mMacrophageLabel = newMacrophageLabel;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::GetCellTypesToEat()
{
    return mCellTypesToEat;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::SetCellTypesToEat(std::vector<unsigned> newCellTypesToEat)
{
    mCellTypesToEat = newCellTypesToEat;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PhenotypeBasedMacrophagePhagocytosisModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    //TODO add parameters here
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}
//
//// Explicit instantiation
//template class PhenotypeBasedMacrophagePhagocytosisModifier<1>;
//template class PhenotypeBasedMacrophagePhagocytosisModifier<2>;
//template class PhenotypeBasedMacrophagePhagocytosisModifier<3>;

// Explicit instantiation
template class PhenotypeBasedMacrophagePhagocytosisModifier<1,1>;
template class PhenotypeBasedMacrophagePhagocytosisModifier<1,2>;
template class PhenotypeBasedMacrophagePhagocytosisModifier<2,2>;
template class PhenotypeBasedMacrophagePhagocytosisModifier<1,3>;
template class PhenotypeBasedMacrophagePhagocytosisModifier<2,3>;
template class PhenotypeBasedMacrophagePhagocytosisModifier<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PhenotypeBasedMacrophagePhagocytosisModifier)

