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

#include "EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs.hpp"
#include "SimpleLinearEllipticSolver.hpp"
//#include "MacrophageCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
		boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
		bool isNeumannBoundaryCondition,
		boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
		double stepSize,
		Vec solution)
		: AbstractBoxDomainPdeModifier<DIM>(pPde,
				pBoundaryCondition,
				isNeumannBoundaryCondition,
				pMeshCuboid,
				stepSize,
				solution)
				{
				}

template<unsigned DIM>
EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::~EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs()
{
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	if ((SimulationTime::Instance()->GetTimeStepsElapsed()%mTimestepInterval)==0)
	{
		// Set up boundary conditions
		std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

		this->UpdateCellPdeElementMap(rCellPopulation);

		// When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
		// Pass in already updated CellPdeElementMap to speed up finding cells.
		this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

		// Use SimpleLinearEllipticSolver as Averaged Source PDE
		///\todo allow other PDE classes to be used with this modifier
		SimpleLinearEllipticSolver<DIM,DIM> solver(this->mpFeMesh,
				boost::static_pointer_cast<AbstractLinearEllipticPde<DIM,DIM> >(this->GetPde()).get(),
				p_bcc.get());

		///\todo Use solution at previous time step as an initial guess for Solve()
		Vec old_solution_copy = this->mSolution;
		this->mSolution = solver.Solve();
		if (old_solution_copy != nullptr)
		{
			PetscTools::Destroy(old_solution_copy);
		}
	}
	this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

	// Call these  methods to solve the PDE on the initial step and output the results
	UpdateAtEndOfTimeStep(rCellPopulation);
	this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

	// To be well-defined, elliptic PDE problems on box domains require at least some Dirichlet boundary conditions
	///\todo Replace this assertion with an exception in the constructor
	assert(!(this->IsNeumannBoundaryCondition()));

	if (!this->mSetBcsOnBoxBoundary) {

		double height = this->mpMeshCuboid->GetWidth(0);
		double width = this->mpMeshCuboid->GetWidth(1);
		if (DIM > 2){
			double depth = this->mpMeshCuboid->GetWidth(2);
		}

		// Loop over the PDE mesh, set dirichlet boundaries on all edges
		for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
			 node_iter != this->mpFeMesh->GetNodeIteratorEnd();
			 ++node_iter) {
			c_vector<double, DIM> pdeNodeLocation = node_iter->rGetLocation();
			bool isVessel = false;
			for (int i = 0; i < mVesselLocations.size(); i++){
				isVessel = isVessel or (pdeNodeLocation[0] == mVesselLocations[i][0] and pdeNodeLocation[1] == mVesselLocations[i][1]);

			};
//			PRINT_VARIABLE(pdeNodeLocation);
//			bool isVessel = (pdeNodeLocation[0] == 5 and pdeNodeLocation[1] == 5)
//						or (pdeNodeLocation[0] == 6 and pdeNodeLocation[1] == 10)
//						or (pdeNodeLocation[0] == 8 and pdeNodeLocation[1] == 21)
//						or (pdeNodeLocation[0] == 15 and pdeNodeLocation[1] == 7)
//						or (pdeNodeLocation[0] == 19 and pdeNodeLocation[1] == 17);
			if (isVessel){
				unsigned node_index = node_iter->GetIndex();
				p_bcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(node_index), this->mpBoundaryCondition.get(),
													 0, false);
			}


//			bool isDirichletBoundaryNode = (0 == pdeNodeLocation[0]);
//			bool isNeumannBoundaryNode = (height == pdeNodeLocation[0] or 0 == pdeNodeLocation[1] or width == pdeNodeLocation[1]);
//			if (isDirichletBoundaryNode){
//				unsigned node_index = node_iter->GetIndex();
//				p_bcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(node_index), this->mpBoundaryCondition.get(),
//													 0, false);
//			}
//			else if (isNeumannBoundaryNode){
//				unsigned node_index = node_iter->GetIndex();
//				std::vector<unsigned> test = this->mpFeMesh->GetContainingElementIndices(pdeNodeLocation);
//				PRINT_VARIABLE(test.size());
//				// Neumann bonudary conditions are set on the element, not the node
//
//				//p_bcc->AddNeumannBoundaryCondition(this->mpFeMesh->GetNode(node_index), this->mpBoundaryCondition.get(),
//				//									 0);
//			}

		}


	}


	return p_bcc;
}

template<unsigned DIM>
int EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::GetTimestepInterval()
{
	return mTimestepInterval;
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::SetTimestepInterval(int newTimestepInterval)
{
	assert(newTimestepInterval >= 0);
	mTimestepInterval = newTimestepInterval;
}


template<unsigned DIM>
std::vector<ChastePoint<DIM> > EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::GetVesselLocations()
{
	return mVesselLocations;
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::SetVesselLocations(std::vector<ChastePoint<DIM> > newVesselLocations)
{
	mVesselLocations = newVesselLocations;
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<TimestepInterval>" << mTimestepInterval << "</TimestepInterval>\n";
	*rParamsFile << "\t\t\t<mVesselLocations>";
	for (int i = 0; i < mVesselLocations.size(); i++){
		*rParamsFile << "[";
		for (int j = 0; j < DIM; j++){
			*rParamsFile << mVesselLocations[i][j] << ",";
		}
		*rParamsFile << "],";
	}
	*rParamsFile << "\t\t\t</mVesselLocations>\n";

	// No parameters to output, so just call method on direct parent class
	AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}



// Explicit instantiation
template class EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<1>;
template class EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<2>;
template class EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs)

