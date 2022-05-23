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

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "ExecutableSupport.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"


#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"


// Write VTU files
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "BoundaryNodeWriter.hpp"

#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"


#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde_edited_phenotypeDependent.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde_edited.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs.hpp"
#include "ApoptoticCellKiller.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellVolumesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellCycleClockWriter.hpp"

#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages.hpp"

#include "CellPropertyCollection.hpp"
#include "CellLabel.hpp"

#include "BoundaryCellWriter.hpp"

#include "GreenspanAndVolumeCellCycle.hpp"
#include "VolumeTrackingModifier.hpp"

#include "DiffusionForceChooseD.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "Arwert2018_MacrophageCellCycle.hpp"
#include "Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier.hpp"
#include "PhenotypeBasedMacrophagePhagocytosisModifier.hpp"
#include "ChemotacticForce_SpecifyNutrientAndCellType.hpp"
#include "ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>


/*
 * Prototype functions
 */
void SetupSingletons(unsigned seed);
void DestroySingletons();
void SetupAndRunSimulation(std::string id_string, double avgCellCycleDurationStroma, double avgCellCycleDurationTumour, double macrophageToCSF1_sensitivity, double macrophageToCXCL12_sensitivity, double tumourCellToEGF_sensitivity, double macrophagePhenotypeIncrementPerHour, double tgfThresholdForPhenotypeSwitch, double distanceToBloodVessels, int numberOfBloodVessels, double halfMaximalExtravasationCsf1Conc, double maximumProbOfExtravasationPerHour, double iterationNumber);
void OutputToConsole(std::string id_string, std::string leading);

int main(int argc, char *argv[])
{
	// This sets up PETSc and prints out copyright information, etc.
	//ExecutableSupport::StandardStartup(&argc, &argv);
	ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

	// Define command line options
	boost::program_options::options_description general_options("This is a Chaste executable.\n");
	general_options.add_options()
            ("help", "produce help message")
            ("ID", boost::program_options::value<std::string>(),"ID string for the simulation")
            ("ACCD_T", boost::program_options::value<double>()->default_value(24),"averageCellCycleDuration_Tumour: Average cell cycle duration for tumour cells")
			("ACCD_S", boost::program_options::value<double>()->default_value(32),"averageCellCycleDuration_Stroma: Average cell cycle duration for stromal cells")
            ("CS_M2CSF", boost::program_options::value<double>()->default_value(2.5),"chemotacticSensitivity_macrophageToCSF: chi value for macrophage chemotaxis sensitivity to CSF")
			("CS_M2CXCL12", boost::program_options::value<double>()->default_value(10),"chemotacticSensitivity_macrophageToCXCL12: chi value for macrophage chemotaxis sensitivity to CXCL12")
			("CS_TC2EGF", boost::program_options::value<double>()->default_value(10),"chemotacticSensitivity_tumourCellToEGF: chi value for tumour cell chemotaxis sensitivity to EGF")
			("MPI", boost::program_options::value<double>()->default_value(2.5),"macrophagePhenotypeIncrementPerHour: amount to increment macrophage phenotype by each hour when above TGF threshold")
			("DTBV", boost::program_options::value<double>()->default_value(7.0),"distanceToBloodVessels: initial distance from tumour cells to blood vessels")
			("NOBV", boost::program_options::value<int>()->default_value(35),"numberOfBloodVessels: target number of blood vessels to include in the simulation")
			("HMECSF", boost::program_options::value<double>()->default_value(0.5),"halfMaximalExtravasationCsf1Conc: concentration of CSF1 at which the extravasation rate will be half maximal")
			("MPOE", boost::program_options::value<double>()->default_value(0.05),"maximumProbOfExtravasationPerHour: maximum rate of extravasation from a blood vessel (infinite CSF1 limit)")
			("TGF_TFPS", boost::program_options::value<double>()->default_value(2.5),"tgfThresholdForPhenotypeSwitch: TGF threshold above which macrophage will begin to polarise")
			("IT", boost::program_options::value<double>(),"Iteration Number");


	// Define parse command line into variables_map
	boost::program_options::variables_map variables_map;
	boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

	// Print help message if wanted
	if (variables_map.count("help"))
	{
		std::cout << setprecision(3) << general_options << "\n";
		std::cout << general_options << "\n";
		return 1;
	}

	// Get ID and name from command line
	std::string id_string = variables_map["ID"].as<std::string>();

	double avgCellCycleDurationStroma = variables_map["ACCD_S"].as<double>();
	double avgCellCycleDurationTumour = variables_map["ACCD_T"].as<double>();
	double macrophageToCSF1_sensitivity = variables_map["CS_M2CSF"].as<double>();
	double macrophageToCXCL12_sensitivity = variables_map["CS_M2CXCL12"].as<double>();
	double tumourCellToEGF_sensitivity = variables_map["CS_TC2EGF"].as<double>();
	double macrophagePhenotypeIncrementPerHour = variables_map["MPI"].as<double>();
	double tgfThresholdForPhenotypeSwitch = variables_map["TGF_TFPS"].as<double>();
	double distanceToBloodVessels = variables_map["DTBV"].as<double>();
	int numberOfBloodVessels = variables_map["NOBV"].as<int>();
	double halfMaximalExtravasationCsf1Conc = variables_map["HMECSF"].as<double>();
	double maximumProbOfExtravasationPerHour = variables_map["MPOE"].as<double>();
	double iterationNumber = variables_map["IT"].as<double>();


	OutputToConsole(id_string, "Started");

	// AWS simulations are seeded randomly
	unsigned seed = time(NULL);

	SetupSingletons(seed);
	SetupAndRunSimulation(id_string, avgCellCycleDurationStroma, avgCellCycleDurationTumour, macrophageToCSF1_sensitivity, macrophageToCXCL12_sensitivity, tumourCellToEGF_sensitivity, macrophagePhenotypeIncrementPerHour, tgfThresholdForPhenotypeSwitch, distanceToBloodVessels, numberOfBloodVessels, halfMaximalExtravasationCsf1Conc, maximumProbOfExtravasationPerHour, iterationNumber);
	DestroySingletons();

	OutputToConsole(id_string, "Completed");
}

void SetupSingletons(unsigned seed)
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);
	//RandomNumberGenerator::Instance()->Reseed(time(NULL));
	//std::stringstream message;
	//message << "Reseeding with seed " << std::to_string(seed) << std::endl;
	//std::cout << message.str() << std::flush;
	RandomNumberGenerator::Instance()->Reseed(seed);
	CellPropertyRegistry::Instance()->Clear();
	CellId::ResetMaxCellId();
}

void DestroySingletons()
{
	// This is from the tearDown method of the test suite
	SimulationTime::Destroy();
	RandomNumberGenerator::Destroy();
	CellPropertyRegistry::Instance()->Clear();
}

void OutputToConsole(std::string id_string, std::string leading)
{
	// Compose the message
	std::stringstream message;
	message << leading << " simulation with ID string " << id_string << std::endl;

	// Send it to the console
	std::cout << message.str() << std::flush;
}

void SetupAndRunSimulation(std::string id_string, double avgCellCycleDurationStroma, double avgCellCycleDurationTumour, double macrophageToCSF1_sensitivity, double macrophageToCXCL12_sensitivity, double tumourCellToEGF_sensitivity, double macrophagePhenotypeIncrementPerHour, double tgfThresholdForPhenotypeSwitch, double distanceToBloodVessels, int numberOfBloodVessels, double halfMaximalExtravasationCsf1Conc, double maximumProbOfExtravasationPerHour, double iterationNumber) {
	{

		int visualisationOutputFrequencyPerHour = 2;
		const int DIM = 2;

		// Generate Mesh:
		// Make Vector
		std::vector<Node<DIM>*> nodes;
		// Add some nodes

		// Domain size
		double domainHeight = 50;
		double domainWidth = 50;

		// Fill the domain with stromal cells

		// Stroma Nodes
		unsigned nodeNum=0;
		unsigned numStromaNodes = 0;
		// For staggering rows
		unsigned rownum = 1;
		double offset;
		for (double x=0; x<=domainWidth; x=x+0.75)//x++)
		{
			offset = 0;
			if (rownum % 2)
			{
				offset = 0.75/2;
			}
			for (double y=0+offset; y<=domainHeight; y=y+0.75)//y++)
			{
				if (DIM == 2) {
					nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
					nodeNum++;
					numStromaNodes++;
				}
			}
			rownum++;
		}
		// Tumour Nodes
		unsigned numTumourNodes = 0;
		// Seed a cluster of tumour cells in the centre
		double x = domainWidth*0.5;
		double y = domainHeight*0.5;
		nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
		nodeNum++;
		nodes.push_back(new Node<DIM>(nodeNum, false, x-0.5, y));
		nodeNum++;
		nodes.push_back(new Node<DIM>(nodeNum, false, x, y-0.5));
		nodeNum++;
		nodes.push_back(new Node<DIM>(nodeNum, false, x-0.5, y-0.5));
		nodeNum++;
		numTumourNodes = numTumourNodes+4;


		// Make cell pointers
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

		unsigned cancerLabelColour = 3;
		unsigned macrophageLabelColour = 7;
		MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
		MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
		for (unsigned i=0; i<numStromaNodes; i++)
		{
			GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
			p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
			p_model->SetHypoxicConcentration(0.1);
			p_model->SetNecroticConcentration(0.01);
			p_model->SetQuiescenceVolumeProportion(0.75);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->GetCellData()->SetItem("csf1", 0);
			p_cell->GetCellData()->SetItem("egf", 0);
			p_cell->GetCellData()->SetItem("cxcl12", 0);
			p_cell->GetCellData()->SetItem("tgf", 0);
			p_cell->GetCellData()->SetItem("phenotype", -2);
			p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
			p_cell->AddCellProperty(p_stroma_label);
			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
							   (avgCellCycleDurationStroma*0.75);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);

		}

		for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
		{

			GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
			p_model->SetDimension(DIM);
			p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
			p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
			p_model->SetHypoxicConcentration(0.01);
			p_model->SetNecroticConcentration(0.01);
			// Tumour cells have much higher tolerance to contact inhibition
			p_model->SetQuiescenceVolumeProportion(0.6);

			CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
			p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("oxygen", 1);
			p_cell->GetCellData()->SetItem("csf1", 0);
			p_cell->GetCellData()->SetItem("egf", 0);
			p_cell->GetCellData()->SetItem("cxcl12", 0);
			p_cell->GetCellData()->SetItem("tgf", 0);
			p_cell->GetCellData()->SetItem("phenotype", -1);
			p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
			p_cell->AddCellProperty(p_cancer_label);

			double birthTime = - RandomNumberGenerator::Instance()->ranf() *
							   (avgCellCycleDurationTumour*0.75);
			p_cell->SetBirthTime(birthTime);
			cells.push_back(p_cell);
		}



		NodesOnlyMesh<DIM> mesh;
		// Cut off length: 1.5 cell radii
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);


		// Make cell population (2D)
		NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

		// Write summary statistic files
		cell_population.AddPopulationWriter<NodeLocationWriter>();
		cell_population.AddPopulationWriter<BoundaryNodeWriter>();
		cell_population.AddCellWriter<CellProliferativeTypesWriter>();
		cell_population.AddCellWriter<BoundaryCellWriter>();
		cell_population.AddCellWriter<CellCycleClockWriter>();


		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		ChastePoint<DIM> lower(0, 0, 0);
		ChastePoint<DIM> upper(domainWidth,domainHeight, 0);
		MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


		//
		// OXYGEN
		//

		// We have Neumann BCs along the boundaries of the simulation domain
		// Dirichlet BCs are set at blood vessel point sources
		// Choose BVs to be random mesh points, rather than picking random points and then meshing
		std::vector<ChastePoint<DIM> > vesselLocations;
        std::vector<ChastePoint<DIM> > possibleVesselList;
		// Choose blood vessel coordinates

		double centreX = 0.5*domainWidth;
		double centreY = 0.5*domainHeight;

		for (int x = 1; x < domainWidth; x++){
			for (int y = 1; y < domainHeight; y++){
				double dist = sqrt(pow((x-centreX),2) + pow((y - centreY),2));
                if ( dist >= distanceToBloodVessels )
                {
                    possibleVesselList.push_back(ChastePoint<DIM> (x,y));
                }
			}
		}

		std::vector<unsigned> order;
		RandomNumberGenerator::Instance()->Shuffle(possibleVesselList.size(), order);
		// Now select numberOfBloodVessels points randomly from possibleVesselList
		// If possibleVesselList contains fewer than numberOfBloodVessels elements, select them all
		if (possibleVesselList.size() < numberOfBloodVessels)
		{
			numberOfBloodVessels = possibleVesselList.size();
		}
		for (int i = 0; i < numberOfBloodVessels; i++){
			vesselLocations.push_back(possibleVesselList[order[i]]);
		}




		double consumptionRate = -0.03;//was -0.03
		double diffusionCoefficient = 1;
		MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
		bool is_neumann_bc = false; // Dirichlet BCs
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		int updateIntervalForPdeInTimesteps = 120/2;
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
		p_oxygen_pde_modifier->SetTimestepInterval(1);
		p_oxygen_pde_modifier->SetVesselLocations(vesselLocations);
		p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
		p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false


		//
		// CXCL12
		//
		double cxcl12_decay_rate = -0.02;
		MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
		bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
		// Create a PDE modifier and set the name of the dependent variable in the PDE
		MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
		p_cxcl12_pde_modifier->SetTimestepInterval(1);
		p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
		p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
		p_cxcl12_pde_modifier->SetOutputGradient(true);
		p_cxcl12_pde_modifier->SetVesselLocations(vesselLocations);
		//
		// CSF-1
		//
		double dudtCoefficient = 1.0;
		double csfDiffusionCoefficient = 1.0;
		double sourceCoefficient = 0.25;
		double decayRate_csf = 0.02;
		MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
		p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
		p_csf1_pde->SetDecayRate(decayRate_csf);
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
		bool is_neumann_bc_csf1 = true;
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
		p_csf1_pde_modifier->SetDependentVariableName("csf1");
		p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
		p_csf1_pde_modifier->SetOutputGradient(true);

		//
		// TGF
		//
		double dudtCoefficient_tgf = 1.0;
		double diffusionCoefficient_tgf = 0.1;
		double sourceCoefficient_tgf = 0.1;
		double decayRate_tgf = 0.1;
		MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
		p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
		p_tgf_pde->SetDecayRate(decayRate_tgf);
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
		bool is_neumann_bc_tgf = true;
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
		p_tgf_pde_modifier->SetDependentVariableName("tgf");
		p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
		p_tgf_pde_modifier->SetOutputGradient(true);

		//
		// EGF
		//
		double dudtCoefficient_egf = 1.0;
		double diffusionCoefficient_egf = 0.2;
		double sourceCoefficient_egf = 0.2;
		double decayRate_egf = 0.1;
		MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited_phenotypeDependent<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
		p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
		p_egf_pde->SetDecayRate(decayRate_egf);
		MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
		bool is_neumann_bc_egf = true;
		MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 1.0));
		p_egf_pde_modifier->SetDependentVariableName("egf");
		p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
		p_egf_pde_modifier->SetOutputGradient(true);



		/* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
		OffLatticeSimulation<DIM> simulator(cell_population);
		simulator.AddSimulationModifier(p_oxygen_pde_modifier);
		simulator.AddSimulationModifier(p_csf1_pde_modifier);
		simulator.AddSimulationModifier(p_egf_pde_modifier);
		simulator.AddSimulationModifier(p_tgf_pde_modifier);
		simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
		std::stringstream output_directory;
		output_directory << "Arwert2018_AWS/" << id_string;
		simulator.SetOutputDirectory(output_directory.str());
		simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
		double burnInPeriod = 0;
		double infiltrationDuration = 500;
		double endTime = burnInPeriod+infiltrationDuration;
		simulator.SetEndTime(endTime);


		MAKE_PTR(Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>, p_addMacs_modifier);
		p_addMacs_modifier->SetTgfbThresholdForMacrophagePhenotypeSwitch(tgfThresholdForPhenotypeSwitch);
		p_addMacs_modifier->SetHalfMaximalExtravasationCsf1Conc(halfMaximalExtravasationCsf1Conc);
		p_addMacs_modifier->SetMaximalProbOfExtravasationPerHour(maximumProbOfExtravasationPerHour);
		p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(infiltrationDuration);
		p_addMacs_modifier->SetVesselLocations(vesselLocations);
		p_addMacs_modifier->SetMacrophagePhenotypeIncrementPerHour(macrophagePhenotypeIncrementPerHour);
		simulator.AddSimulationModifier(p_addMacs_modifier);

        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
        simulator.AddSimulationModifier(p_phagoytosis_modifier);


		// Add Brownian motion for macrophages
		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceMac);
		p_diffusion_forceMac->SetDiffusionScalingConstant(0.01);
		p_diffusion_forceMac->SetCellTypesForDiffusion({macrophageLabelColour});
		simulator.AddForce(p_diffusion_forceMac);

		MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
		simulator.AddSimulationModifier(p_modifier);

		// Macrophages go up csf1 gradient
		MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
		p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
		p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
		p_chemotaxis_MacrophageToCSF1->SetSensitivity(macrophageToCSF1_sensitivity);
		simulator.AddForce(p_chemotaxis_MacrophageToCSF1);

		// Tumour cells go up egf gradient
		MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
		p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
		p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
		p_chemotaxis_TumourCellToEGF->SetSensitivity(tumourCellToEGF_sensitivity);
		simulator.AddForce(p_chemotaxis_TumourCellToEGF);

		// Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
		MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
		p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
		p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
		p_chemotaxis_MacrophageToCXCL12->SetSensitivity(macrophageToCXCL12_sensitivity);
		simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


		MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
		p_force->SetMeinekeSpringStiffness(5.0);
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);

		// Killer which removes apoptotic cells
		MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
		simulator.AddCellKiller(p_apoptosis_killer);

		// Boundary conditions: Add a wall along the x axis
		c_vector<double,2> point = zero_vector<double>(2);
		c_vector<double,2> normal = zero_vector<double>(2);
		normal(0) = -1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);

		// Boundary conditions: Add a wall along the bottom
		c_vector<double,2> normal2 = zero_vector<double>(2);
		c_vector<double,2> point2 = zero_vector<double>(2);
		normal2(1) = -1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
		simulator.AddCellPopulationBoundaryCondition(p_bc2);

		// Boundary conditions: Add a wall along the top
		c_vector<double,2> normal3 = zero_vector<double>(2);
		c_vector<double,2> point3 = zero_vector<double>(2);
		normal3(1) = 1.0;
		point3(1) = domainHeight;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
		simulator.AddCellPopulationBoundaryCondition(p_bc3);

		// Boundary conditions: Add a wall along the right hand side
		c_vector<double,2> point4 = zero_vector<double>(2);
		c_vector<double,2> normal4 = zero_vector<double>(2);
		point4(0) = domainWidth;
		normal4(0) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
		simulator.AddCellPopulationBoundaryCondition(p_bc4);



		/* To run the simulation, we call {{{Solve()}}}. */
		simulator.Solve();

		/* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}

	}
}
