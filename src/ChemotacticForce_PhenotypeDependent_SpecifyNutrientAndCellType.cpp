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

#include "ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType.hpp"

#include "CellLabel.hpp"

template<unsigned DIM>
ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType()
    : AbstractForce<DIM>(),
            mNutrient("nutrient"),
            mCellTypesForChemotaxis({0}),
            mSensitivity(1)
{
}

template<unsigned DIM>
ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::~ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType()
{
}

template<unsigned DIM>
double ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::GetChemotacticForceMagnitude(const double concentration, const double concentrationGradientMagnitude)
{
    return mSensitivity;//concentration; // temporary force law - can be changed to something realistic
                          // without tests failing
}

template<unsigned DIM>
void ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only labelled cells move chemotactically
        CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection();
        unsigned label = 0;
        if (collection.HasProperty<CellLabel>())
        {
            CellPropertyCollection collectionLabel = collection.GetProperties<CellLabel>();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collectionLabel.GetProperty());
            label = p_label->GetColour();
        }

        bool isCellSubjectToChemotaxis = std::find(mCellTypesForChemotaxis.begin(), mCellTypesForChemotaxis.end(), label) != mCellTypesForChemotaxis.end();
        if (isCellSubjectToChemotaxis)
        {
            unsigned node_global_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double,DIM> gradient;
            // Find gradients in x, y, z directions at node
            if (DIM == 2){
                gradient(0) = cell_iter->GetCellData()->GetItem(mNutrient + "_grad_x");
                gradient(1) = cell_iter->GetCellData()->GetItem(mNutrient + "_grad_y");
            }
            else if (DIM == 3){
                gradient(0) = cell_iter->GetCellData()->GetItem(mNutrient + "_grad_x");
                gradient(1) = cell_iter->GetCellData()->GetItem(mNutrient + "_grad_y");
                gradient(2) = cell_iter->GetCellData()->GetItem(mNutrient + "_grad_z");
            }

            // Need to calculate the gradient
            double nutrient_concentration = cell_iter->GetCellData()->GetItem(mNutrient);
            double magnitude_of_gradient = norm_2(gradient);

            double force_magnitude = GetChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);
            double phenotype = cell_iter->GetCellData()->GetItem("phenotype");

            //todo replace with something more refined
            if (mNutrient == "cxcl12") {
                force_magnitude = force_magnitude * phenotype;
            }
            else if (mNutrient == "csf1") {
                force_magnitude = force_magnitude * (1-phenotype);
            }
            else
                {
                assert(1==2);
            }
            // force += chi * gradC/|gradC|
            if (magnitude_of_gradient > 0)
            {
                c_vector<double,DIM> force = force_magnitude*(gradient/magnitude_of_gradient);//(force_magnitude/magnitude_of_gradient)*gradient;
                rCellPopulation.GetNode(node_global_index)->AddAppliedForceContribution(force);
            }
            // else Fc=0


        }
    }
}

template<unsigned DIM>
void ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template <unsigned DIM>
std::vector<unsigned> ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::GetCellTypesForChemotaxis()
{
    return mCellTypesForChemotaxis;
}

template <unsigned DIM>
void ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::SetCellTypesForChemotaxis(std::vector<unsigned> newCellTypesForChemotaxis)
{
    mCellTypesForChemotaxis = newCellTypesForChemotaxis;
}

template <unsigned DIM>
std::string ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::GetNutrient()
{
    return mNutrient;
}

template <unsigned DIM>
void ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::SetNutrient(std::string newNutrient)
{
    mNutrient = newNutrient;
}

template <unsigned DIM>
double ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::GetSensitivity()
{
    return mSensitivity;
}

template <unsigned DIM>
void ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>::SetSensitivity(double newSensitivity)
{
    mSensitivity = newSensitivity;
}

// Explicit instantiation
template class ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<1>;
template class ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<2>;
template class ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType)
