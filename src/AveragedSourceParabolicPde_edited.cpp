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

#include "AveragedSourceParabolicPde_edited.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
AveragedSourceParabolicPde_edited<DIM>::AveragedSourceParabolicPde_edited(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double duDtCoefficient,
                                                            double diffusionCoefficient,
                                                            double sourceCoefficient)
    : AveragedSourceParabolicPde<DIM>(rCellPopulation,duDtCoefficient,diffusionCoefficient,sourceCoefficient),
      mCellTypesAsSource({0}),
      mDecayRate(0.1)
//    mrCellPopulation(rCellPopulation),
//      mDuDtCoefficient(duDtCoefficient),
//      mDiffusionCoefficient(diffusionCoefficient),
//      mSourceCoefficient(sourceCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& AveragedSourceParabolicPde_edited<DIM>::rGetCellPopulation() const
{
    return this->mrCellPopulation;
}

template<unsigned DIM>
void AveragedSourceParabolicPde_edited<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    this->mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<this->mCellDensityOnCoarseElements.size(); elem_index++)
    {
        this->mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned elem_index = 0;
        const ChastePoint<DIM>& r_position_of_cell = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

        if (pCellPdeElementMap != nullptr)
        {
            elem_index = (*pCellPdeElementMap)[*cell_iter];
        }
        else
        {
            elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
        }

        // Update element map if cell has moved
        bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();

        CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection();
        unsigned label = 0;
        if (collection.HasProperty<CellLabel>())
        {
            CellPropertyCollection collectionLabel = collection.GetProperties<CellLabel>();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collectionLabel.GetProperty());
            label = p_label->GetColour();
        }

        bool isTargetCellASource = std::find(mCellTypesAsSource.begin(), mCellTypesAsSource.end(), label) != mCellTypesAsSource.end();
        if (isTargetCellASource and !cell_is_apoptotic)
        {
            double increment = 1.0;
            // If this is phenotype dependent, instead of incrementing by 1.0 we have a step function based on macrophage phenotype

            this->mCellDensityOnCoarseElements[elem_index] += increment;
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<this->mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        this->mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

template<unsigned DIM>
double AveragedSourceParabolicPde_edited<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return this->mDuDtCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde_edited<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    assert(!this->mCellDensityOnCoarseElements.empty());
    double coefficient = this->mSourceCoefficient * this->mCellDensityOnCoarseElements[pElement->GetIndex()];
    // The source term is rho(x)*mSourceCoefficient - mDecayRate*u
    return coefficient - mDecayRate*u;
}

// LCOV_EXCL_START
template<unsigned DIM>
double AveragedSourceParabolicPde_edited<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
c_matrix<double,DIM,DIM> AveragedSourceParabolicPde_edited<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return this->mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double AveragedSourceParabolicPde_edited<DIM>::GetUptakeRateForElement(unsigned elementIndex)
{
    this->mCellDensityOnCoarseElements[elementIndex];
    return this->mCellDensityOnCoarseElements[elementIndex];
}

template <unsigned DIM>
std::vector<unsigned> AveragedSourceParabolicPde_edited<DIM>::GetCellTypesAsSource()
{
    return mCellTypesAsSource;
}

template <unsigned DIM>
void AveragedSourceParabolicPde_edited<DIM>::SetCellTypesAsSource(std::vector<unsigned> newCellTypesAsSource)
{
    mCellTypesAsSource = newCellTypesAsSource;
}

template <unsigned DIM>
double AveragedSourceParabolicPde_edited<DIM>::GetDecayRate()
{
    return mDecayRate;
}

template <unsigned DIM>
void AveragedSourceParabolicPde_edited<DIM>::SetDecayRate(double newDecayRate)
{
    mDecayRate = newDecayRate;
}

// Explicit instantiation
template class AveragedSourceParabolicPde_edited<1>;
template class AveragedSourceParabolicPde_edited<2>;
template class AveragedSourceParabolicPde_edited<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AveragedSourceParabolicPde_edited)
