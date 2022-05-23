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

#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages.hpp"
#include "CellLabel.hpp"
#include "MacrophageCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages()
: GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
  mHomotypicLabelledSpringConstantMultiplier(1.0),
  mHeterotypicLabelledSpringConstantMultiplier(1.0),
  mCloserThanRestLengthSpringConstantMultiplier(1.0),
  mMacrophageContactRadius(1.5)
  {
  }

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::AreBothCellsInContactWithAMacrophage(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{

	auto nbcp = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>* >(&rCellPopulation);

	std::set<unsigned> neighbours_NodeA = nbcp->GetNodesWithinNeighbourhoodRadius(nodeAGlobalIndex, mMacrophageContactRadius);
	std::set<unsigned> neighbours_NodeB = nbcp->GetNodesWithinNeighbourhoodRadius(nodeBGlobalIndex, mMacrophageContactRadius);

	std::vector<unsigned> sharedNeighbours;
	std::set_intersection(neighbours_NodeA.begin(),neighbours_NodeA.end(),neighbours_NodeB.begin(),neighbours_NodeB.end(),std::back_inserter(sharedNeighbours));

	if (sharedNeighbours.size() > 0)
	{

		// Now check each index to see if a shared neighbour is a macrophage
		for(unsigned i = 0;i<sharedNeighbours.size();i++)
		{

			unsigned candidateIndex = sharedNeighbours[i];
			CellPtr p_sharedCell = rCellPopulation.GetCellUsingLocationIndex(candidateIndex);

			if (p_sharedCell->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>())
			{

				return true;
			}
		}
	}

	return false;

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
		unsigned nodeAGlobalIndex,
		unsigned nodeBGlobalIndex,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
		bool isCloserThanRestLength)
		{

	if (isCloserThanRestLength)
	{
		return mCloserThanRestLengthSpringConstantMultiplier;
	}
	else
	{
		double a_multiplier;
		double b_multiplier;
		double time_until_death_a;
		double time_until_death_b;

		// Determine which (if any) of the cells corresponding to these nodes are labelled apoptotic
		CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
		bool cell_A_is_labelled = p_cell_A->HasApoptosisBegun();

		CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);
		bool cell_B_is_labelled = p_cell_B->HasApoptosisBegun();

		if(cell_A_is_labelled)
		{
			time_until_death_a = p_cell_A->GetTimeUntilDeath();
			a_multiplier = time_until_death_a / p_cell_A->GetApoptosisTime();

			if(cell_B_is_labelled)
			{
				time_until_death_b = p_cell_B->GetTimeUntilDeath();
				b_multiplier = time_until_death_b / p_cell_B->GetApoptosisTime();

				return a_multiplier*b_multiplier*mHomotypicLabelledSpringConstantMultiplier;
			}
			else
			{
				return a_multiplier*mHeterotypicLabelledSpringConstantMultiplier;
			}
		}
		else
		{
			if(cell_B_is_labelled)
			{
				time_until_death_b = p_cell_B->GetTimeUntilDeath();
				b_multiplier = time_until_death_b / p_cell_B->GetApoptosisTime();

				return b_multiplier*mHeterotypicLabelledSpringConstantMultiplier;
			}
			else
			{
				// If neither cell is apoptotic, check for microbead interactions
				if((p_cell_A->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>()) and (p_cell_B->GetCellProliferativeType()->template IsType<MacrophageCellProliferativeType>()))
				{
					return 0.0;
				}
				else
				{

					// No bead-bead means that spring remains unchanged
					if(AreBothCellsInContactWithAMacrophage(nodeAGlobalIndex,nodeBGlobalIndex,rCellPopulation))
					{

						return 0; // Arbitrary, change this
					}
					else
					{

						return 1.0;
					}
				}

			}
		}
	}
		}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::GetMacrophageContactRadius()
{
	return mMacrophageContactRadius;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::SetMacrophageContactRadius(double newMacrophageContactRadius)
{
	mMacrophageContactRadius = newMacrophageContactRadius;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::GetHomotypicLabelledSpringConstantMultiplier()
{
	return mHomotypicLabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::SetHomotypicLabelledSpringConstantMultiplier(double labelledSpringConstantMultiplier)
{
	assert(labelledSpringConstantMultiplier > 0.0);
	mHomotypicLabelledSpringConstantMultiplier = labelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::GetHeterotypicLabelledSpringConstantMultiplier()
{
	return mHeterotypicLabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::SetHeterotypicLabelledSpringConstantMultiplier(double heterotypicLabelledSpringConstantMultiplier)
{
	assert(heterotypicLabelledSpringConstantMultiplier > 0.0);
	mHeterotypicLabelledSpringConstantMultiplier = heterotypicLabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::GetCloserThanRestLengthSpringConstantMultiplier()
{
	return mCloserThanRestLengthSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::SetCloserThanRestLengthSpringConstantMultiplier(double closerThanRestLengthSpringConstantMultiplier)
{
	assert(closerThanRestLengthSpringConstantMultiplier > 0.0);
	mCloserThanRestLengthSpringConstantMultiplier = closerThanRestLengthSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<HomotypicLabelledSpringConstantMultiplier>" << mHomotypicLabelledSpringConstantMultiplier << "</HomotypicLabelledSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<mHeterotypicLabelledSpringConstantMultiplier>" << mHeterotypicLabelledSpringConstantMultiplier << "</mHeterotypicLabelledSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<mCloserThanRestLengthSpringConstantMultiplier>" << mCloserThanRestLengthSpringConstantMultiplier << "</mCloserThanRestLengthSpringConstantMultiplier>\n";

	// Call direct parent class
	GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<1,1>;
template class GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<1,2>;
template class GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<2,2>;
template class GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<1,3>;
template class GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<2,3>;
template class GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages)
