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

#ifndef ARWERT2018_ADDMACROPHAGESFROMPOINTVESSELINRESPONSETOCSFMODIFIER_HPP
#define ARWERT2018_ADDMACROPHAGESFROMPOINTVESSELINRESPONSETOCSFMODIFIER_HPP

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"


template<unsigned DIM>
class Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

protected:

    double timeToAddMacrophages;

    unsigned numberOfMacrophagesToAdd;

    bool haveMacrophagesBeenAdded;

    unsigned numberOfHoursToRunSimulationAfterAddingMacrophages;

    double macrophagePhenotypeIncrementPerHour;

    double tgfbThresholdForMacrophagePhenotypeSwitch;

    double halfMaximalExtravasationCsf1Conc;

    double maximalProbOfExtravasationPerHour;

    /* Locations of points for blood vessels */
    std::vector<ChastePoint<DIM> > mVesselLocations;

public:

    /**
     * Default constructor.
     */
    Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier();

    /**
     * Destructor.
     */
    virtual ~Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	/**
	 * @return #numberOfTumourCellsBeforeMacrophagesAreAdded.
	 */
    double GetTimeToAddMacrophages();

	/**
	 * Set numberOfTumourCellsBeforeMacrophagesAreAdded.
	 *
	 * @param newNumberOfTumourCellsBeforeMacrophagesAreAdded the new value of timestep
	 */
	void SetTimeToAddMacrophages(double newTimeToAddMacrophages);

	/**
	 * @return #numberOfMacrophagesToAdd.
	 */
	unsigned GetNumberOfMacrophagesToAdd();

	/**
	 * Set numberOfMacrophagesToAdd.
	 *
	 * @param newNumberOfMacrophagesToAdd the new value of timestep
	 */
	void SetNumberOfMacrophagesToAdd(unsigned newNumberOfMacrophagesToAdd);

    /**
     * @return #tgfbThresholdForMacrophagePhenotypeSwitch.
     */
    double GetTgfbThresholdForMacrophagePhenotypeSwitch();

    /**
     * Set tgfbThresholdForMacrophagePhenotypeSwitch.
     *
     * @param newTgfbThresholdForMacrophagePhenotypeSwitch the new value of timestep
     */
    void SetTgfbThresholdForMacrophagePhenotypeSwitch(double newTgfbThresholdForMacrophagePhenotypeSwitch);


    /**
     * @return #macrophagePhenotypeIncrementPerHour.
     */
    double GetMacrophagePhenotypeIncrementPerHour();

    /**
     * Set numberOfMacrophagesToAdd.
     *
     * @param newMacrophagePhenotypeIncrementPerHour the new value of timestep
     */
    void SetMacrophagePhenotypeIncrementPerHour(double newMacrophagePhenotypeIncrementPerHour);

	/**
	 * @return #numberOfHoursToRunSimulationAfterAddingMacrophages.
	 */
	unsigned GetNumberOfHoursToRunSimulationAfterAddingMacrophages();

	/**
	 * Set numberOfHoursToRunSimulationAfterAddingMacrophages.
	 *
	 * @param newNumberOfHoursToRunSimulationAfterAddingMacrophages the new value of timestep
	 */
	void SetNumberOfHoursToRunSimulationAfterAddingMacrophages(unsigned newNumberOfHoursToRunSimulationAfterAddingMacrophages);
    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    std::vector<ChastePoint<DIM> > GetVesselLocations();

    void SetVesselLocations(std::vector<ChastePoint<DIM> > newVesselLocations);

    double GetHalfMaximalExtravasationCsf1Conc();

    void SetHalfMaximalExtravasationCsf1Conc(double newHalfMaximalExtravasationCsf1Conc);

    double GetMaximalProbOfExtravasationPerHour();

    void SetMaximalProbOfExtravasationPerHour(double newMaximalProbOfExtravasationPerHour);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier)

#endif /*ARWERT2018_ADDMACROPHAGESFROMPOINTVESSELINRESPONSETOCSFMODIFIER_HPP*/
