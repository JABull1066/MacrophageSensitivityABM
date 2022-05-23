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

#include "DiffusionForceChooseD.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellLabel.hpp"
//#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
DiffusionForceChooseD<DIM>::DiffusionForceChooseD()
    : AbstractForce<DIM>(),
      mD(0.0),
      mCellTypesForDiffusion({0})
{
	//MARK;
}

template<unsigned DIM>
DiffusionForceChooseD<DIM>::~DiffusionForceChooseD()
{
}

template<unsigned DIM>
void DiffusionForceChooseD<DIM>::SetDiffusionScalingConstant(double newD)
{
	//MARK;
	mD = newD;
}

template<unsigned DIM>
double DiffusionForceChooseD<DIM>::GetDiffusionScalingConstant()
{
	/*
	 * In DiffusionForce, this is given by msBoltzmannConstant*mAbsoluteTemperature/(6.0*mViscosity*M_PI). This method is designed
	 * so that this can be overwritten with an arbitrary constant. The reason for allowing this (rather than defining in terms
	 * of physical constants) is that in the case of cell motion, we want to be able to set a parameter which roughly corresponds
	 * to "cell motility", which isn't defined in terms of these physical constants. This simplifies
	 * the process of having to select arbitrary "temperature" and "viscosity" for what is an active process (cell movement)
	 * that isn't really diffusion at all.
	 */
    return mD;
}

template <unsigned DIM>
std::vector<unsigned> DiffusionForceChooseD<DIM>::GetCellTypesForDiffusion()
{
    return mCellTypesForDiffusion;
}

template <unsigned DIM>
void DiffusionForceChooseD<DIM>::SetCellTypesForDiffusion(std::vector<unsigned> newCellTypesForDiffusion)
{
    mCellTypesForDiffusion = newCellTypesForDiffusion;
}

template<unsigned DIM>
void DiffusionForceChooseD<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    //std::stringstream message;
    //double time = SimulationTime::Instance()->GetTime();
    //message << "Diffusion force " << RandomNumberGenerator::Instance()->ranf() << ", time " << time << std::endl;
    //std::cout << message.str() << std::flush;
    AbstractOffLatticeCellPopulation<DIM>* p_pop = dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation);
    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter) {
        unsigned node_index = node_iter->GetIndex();

        //todo Check that this is the correct index for getting the cell
        CellPropertyCollection collection = (p_pop->GetCellUsingLocationIndex(
                node_index))->rGetCellPropertyCollection();
        unsigned label = 0;
        if (collection.HasProperty<CellLabel>()) {
            CellPropertyCollection collectionLabel = collection.GetProperties<CellLabel>();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collectionLabel.GetProperty());
            label = p_label->GetColour();
        }

        bool isCellSubjectToDiffusion =
                std::find(mCellTypesForDiffusion.begin(), mCellTypesForDiffusion.end(), label) !=
                mCellTypesForDiffusion.end();
        if (isCellSubjectToDiffusion) {

            // Get the index, radius and damping constant of this node
            double node_radius = node_iter->GetRadius();

            // If the node radius is zero, then it has not been set...
            if (node_radius == 0.0) {
                // ...so throw an exception to avoid dividing by zero when we compute diffusion_constant below
                EXCEPTION(
                        "SetRadius() must be called on each Node before calling DiffusionForce::AddForceContribution() to avoid a division by zero error");
            }

            double nu = p_pop->GetDampingConstant(node_index);

            double D = GetDiffusionScalingConstant();
            double diffusion_constant = D / (node_radius * 2);

            c_vector<double, DIM> force_contribution;
            for (unsigned i = 0; i < DIM; i++) {
                /*
                 * The force on this cell is scaled with the timestep such that when it is
                 * used in the discretised equation of motion for the cell, we obtain the
                 * correct formula
                 *
                 * x_new = x_old + sqrt(2*D*dt)*W
                 *
                 * where W is a standard normal random variable.
                 */
                double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

                //force_contribution[i] = (nu*sqrt(2.0*diffusion_constant*dt)/dt)*xi;
                //We want xi to be normally distributed, with variance 2DT
                force_contribution[i] = (nu * sqrt(2.0 * diffusion_constant * dt) / dt) * xi;
            }
            node_iter->AddAppliedForceContribution(force_contribution);
        }
    }
}

template<unsigned DIM>
void DiffusionForceChooseD<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionCoefficient>" << mD << "</DiffusionCoefficient> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DiffusionForceChooseD<1>;
template class DiffusionForceChooseD<2>;
template class DiffusionForceChooseD<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionForceChooseD)
