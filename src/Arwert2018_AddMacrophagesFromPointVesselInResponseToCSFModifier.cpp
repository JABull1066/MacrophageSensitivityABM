/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier.hpp"
#include "SmartPointers.hpp"
#include "Arwert2018_MacrophageCellCycle.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "ChastePoint.hpp"
#include "RandomNumberGenerator.hpp"
#include "NodesOnlyMesh.hpp"

template<unsigned DIM>
Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier()
: AbstractCellBasedSimulationModifier<DIM>(),
  timeToAddMacrophages(0.0),
  numberOfMacrophagesToAdd(100),
  haveMacrophagesBeenAdded(false),
  macrophagePhenotypeIncrementPerHour(0.1),
  tgfbThresholdForMacrophagePhenotypeSwitch(2.0),
  numberOfHoursToRunSimulationAfterAddingMacrophages(UINT_MAX),
  mVesselLocations(),
  maximalProbOfExtravasationPerHour(0.01),
  halfMaximalExtravasationCsf1Conc(1.0)
  {
  }

template<unsigned DIM>
Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::~Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier()
{
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetTimeToAddMacrophages()
{
	return timeToAddMacrophages;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetTimeToAddMacrophages(double newTimeToAddMacrophages)
{
	timeToAddMacrophages = newTimeToAddMacrophages;
}

template<unsigned DIM>
double Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetMaximalProbOfExtravasationPerHour()
{
    return maximalProbOfExtravasationPerHour;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetMaximalProbOfExtravasationPerHour(double newMaximalProbOfExtravasationPerHour)
{
    maximalProbOfExtravasationPerHour = newMaximalProbOfExtravasationPerHour;
}

template<unsigned DIM>
double Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetHalfMaximalExtravasationCsf1Conc()
{
    return halfMaximalExtravasationCsf1Conc;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetHalfMaximalExtravasationCsf1Conc(double newHalfMaximalExtravasationCsf1Conc)
{
    halfMaximalExtravasationCsf1Conc = newHalfMaximalExtravasationCsf1Conc;
}

template<unsigned DIM>
double Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetTgfbThresholdForMacrophagePhenotypeSwitch()
{
    return tgfbThresholdForMacrophagePhenotypeSwitch;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetTgfbThresholdForMacrophagePhenotypeSwitch(double newTgfbThresholdForMacrophagePhenotypeSwitch)
{
    tgfbThresholdForMacrophagePhenotypeSwitch = newTgfbThresholdForMacrophagePhenotypeSwitch;
}

template<unsigned DIM>
double Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetMacrophagePhenotypeIncrementPerHour()
{
    return macrophagePhenotypeIncrementPerHour;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetMacrophagePhenotypeIncrementPerHour(double newMacrophagePhenotypeIncrementPerHour)
{
    macrophagePhenotypeIncrementPerHour = newMacrophagePhenotypeIncrementPerHour;
}

template<unsigned DIM>
unsigned Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetNumberOfMacrophagesToAdd()
{
	return numberOfMacrophagesToAdd;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetNumberOfMacrophagesToAdd(unsigned newNumberOfMacrophagesToAdd)
{
	numberOfMacrophagesToAdd = newNumberOfMacrophagesToAdd;
}

template<unsigned DIM>
unsigned Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetNumberOfHoursToRunSimulationAfterAddingMacrophages()
{
	return numberOfHoursToRunSimulationAfterAddingMacrophages;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetNumberOfHoursToRunSimulationAfterAddingMacrophages(unsigned newNumberOfHoursToRunSimulationAfterAddingMacrophages)
{
	numberOfHoursToRunSimulationAfterAddingMacrophages = newNumberOfHoursToRunSimulationAfterAddingMacrophages;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);
}



template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{


        // Make sure the cell population is updated
        rCellPopulation.Update();

        // Check each blood vessel
		for (int i = 0; i < mVesselLocations.size(); i++){

		    // Check CSF concentration at this vessel
			// This turns out to be annoyingly difficult
			// Instead, use the CSF1 conc at the closest cell
			double csf1Concentration = 0;
			double closestDistance = DBL_MAX;

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                // Is node within 1 unit of this point?
                c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                double cellDistanceFromVessel = 0;
                for (unsigned j=0; i<DIM; i++)
                {
                    cellDistanceFromVessel+=pow((cell_location[j] - mVesselLocations[j][i]),2);
                }
                cellDistanceFromVessel = sqrt(cellDistanceFromVessel);

                if (cellDistanceFromVessel < closestDistance)
                {
                    closestDistance = cellDistanceFromVessel;
                    csf1Concentration = (*cell_iter)->GetCellData()->GetItem("csf1");
                }
            }

            double probOfExtravasationPerHour = maximalProbOfExtravasationPerHour * csf1Concentration / (csf1Concentration + halfMaximalExtravasationCsf1Conc);

			// Turn that into an extravasation probability
            /*
             * We assume a constant time step and that there are an integer number (n = 1/dt)
             * of time steps per hour. We also assume that this method is called every time step
             * and that the probabilities of dying at different times are independent.
             *
             * Let q=mProbabilityOfDeathInAnHour and p="probability of death in a given time step".
             *
             * Probability of not dying in an hour:
             * (1-q) = (1-p)^n = (1-p)^(1/dt).
             *
             * Rearranging for p:
             * p = 1 - (1-q)^dt.
             */
            double extravasation_prob_this_timestep = 1.0 - pow((1.0 - probOfExtravasationPerHour), SimulationTime::Instance()->GetTimeStep());
			// If probable, extravasate a macrophage

			// Optional: Keep track of the total number of macrophages present so it can be capped if wanted

            if (RandomNumberGenerator::Instance()->ranf() < extravasation_prob_this_timestep)
            {
                auto nbcp = dynamic_cast<NodeBasedCellPopulation<DIM>* >(&rCellPopulation);
                // Make a macrophage

                MAKE_PTR(MacrophageCellProliferativeType, p_macrophage_type);
                MAKE_PTR(WildTypeCellMutationState, p_state);
                unsigned macrophageLabelColour = 7;
                MAKE_PTR_ARGS(CellLabel, p_macrophage_label, (macrophageLabelColour));

                double random_x = mVesselLocations[i][0] + 0.25*RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
                double random_y = mVesselLocations[i][1] + 0.25*RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

                // Make Macrophage
                Arwert2018_MacrophageCellCycle *p_model = new Arwert2018_MacrophageCellCycle;
                p_model->SetDimension(DIM);
                p_model->SetPhagocytosisCooldownPeriod(4.0);
                p_model->SetPhenotypeIncrementPerHour(macrophagePhenotypeIncrementPerHour);
                p_model->SetCriticalTGFbetaThreshold(tgfbThresholdForMacrophagePhenotypeSwitch);

                CellPtr pNewCell(new Cell(p_state, p_model));
                pNewCell->SetCellProliferativeType(p_macrophage_type);
                pNewCell->GetCellData()->SetItem("oxygen", 1);
                pNewCell->GetCellData()->SetItem("csf1", 0);
                pNewCell->GetCellData()->SetItem("csf1_grad_x", 0);
                pNewCell->GetCellData()->SetItem("csf1_grad_y", 0);
                pNewCell->GetCellData()->SetItem("egf", 0);
                pNewCell->GetCellData()->SetItem("egf_grad_x", 0);
                pNewCell->GetCellData()->SetItem("egf_grad_y", 0);
                pNewCell->GetCellData()->SetItem("cxcl12", 0);
                pNewCell->GetCellData()->SetItem("cxcl12_grad_x", 0);
                pNewCell->GetCellData()->SetItem("cxcl12_grad_y", 0);
                pNewCell->GetCellData()->SetItem("tgf", 0);
                pNewCell->GetCellData()->SetItem("tgf_grad_x", 0);
                pNewCell->GetCellData()->SetItem("tgf_grad_y", 0);
                pNewCell->GetCellData()->SetItem("phenotype", 0);
                pNewCell->AddCellProperty(p_macrophage_label);

                // 10000 here is a placeholder - it should be overwritten when we add the node to the nbcp...I hope
                Node<DIM>* p_new_node;
                if(DIM == 1)
                {
                    p_new_node = new Node<DIM>(100000,  false,  random_x);
                }
                if(DIM == 2)
                {
                    p_new_node = new Node<DIM>(100000,  false,  random_x, random_y);
                }
                if(DIM == 3)
                {
                    NEVER_REACHED;
                }

                p_new_node->ClearAppliedForce(); // In case velocity is ouptut on the same timestep as the cell has divided
                p_new_node->SetRadius(0.5); //
                //NodesOnlyMesh<3> mesh = nbcp->rGetMesh();
                unsigned new_node_index = nbcp->rGetMesh().AddNode(p_new_node);


                // Update cells vector
                nbcp->rGetCells().push_back(pNewCell);

                // Update mappings between cells and location indices
                nbcp->SetCellUsingLocationIndex(new_node_index, pNewCell);

            }


        }

}


template<unsigned DIM>
std::vector<ChastePoint<DIM> > Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::GetVesselLocations()
{
	return mVesselLocations;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::SetVesselLocations(std::vector<ChastePoint<DIM> > newVesselLocations)
{
	mVesselLocations = newVesselLocations;
}

template<unsigned DIM>
void Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

	*rParamsFile << "\t\t\t<timeToAddMacrophages>" << timeToAddMacrophages << "</timeToAddMacrophages>\n";
	*rParamsFile << "\t\t\t<numberOfMacrophagesToAdd>" << numberOfMacrophagesToAdd << "</numberOfMacrophagesToAdd>\n";
	*rParamsFile << "\t\t\t<numberOfHoursToRunSimulationAfterAddingMacrophages>" << numberOfHoursToRunSimulationAfterAddingMacrophages << "</numberOfHoursToRunSimulationAfterAddingMacrophages>\n";
    *rParamsFile << "\t\t\t<maximalProbOfExtravasationPerHour>" << maximalProbOfExtravasationPerHour << "</maximalProbOfExtravasationPerHour>\n";
    *rParamsFile << "\t\t\t<halfMaximalExtravasationCsf1Conc>" << halfMaximalExtravasationCsf1Conc << "</halfMaximalExtravasationCsf1Conc>\n";
    *rParamsFile << "\t\t\t<macrophagePhenotypeIncrementPerHour>" << macrophagePhenotypeIncrementPerHour << "</macrophagePhenotypeIncrementPerHour>\n";
    *rParamsFile << "\t\t\t<tgfbThresholdForMacrophagePhenotypeSwitch>" << tgfbThresholdForMacrophagePhenotypeSwitch << "</tgfbThresholdForMacrophagePhenotypeSwitch>\n";


	// Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<1>;
template class Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<2>;
template class Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier)

