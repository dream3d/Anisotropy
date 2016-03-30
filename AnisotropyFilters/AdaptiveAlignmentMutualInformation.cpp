/* ============================================================================
* Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
* Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright notice, this
* list of conditions and the following disclaimer in the documentation and/or
* other materials provided with the distribution.
*
* Neither the name of Czech Academy of Sciences, Institute of Physics, nor the names
* of its contributors may be used to endorse or promote products derived from this
* software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The code contained herein was partially funded by the followig grants:
*    Czech Science Foundation (GA CR), project no. GBP108/12/G043
*    Czech Ministry of Education, Youth and Sports (MSMT), project no. LM2015087
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "AdaptiveAlignmentMutualInformation.h"

#include <fstream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/SIMPLibVersion.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersWriter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/Utilities/SIMPLibRandom.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "Anisotropy/AnisotropyConstants.h"

#include "moc_AdaptiveAlignmentMutualInformation.cpp"
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignmentMutualInformation::AdaptiveAlignmentMutualInformation() :
AdaptiveAlignment(),
m_MisorientationTolerance(5.0f),
m_UseGoodVoxels(true),
m_QuatsArrayPath(DREAM3D::Defaults::ImageDataContainerName, DREAM3D::Defaults::CellAttributeMatrixName, DREAM3D::CellData::Quats),
m_CellPhasesArrayPath(DREAM3D::Defaults::ImageDataContainerName, DREAM3D::Defaults::CellAttributeMatrixName, DREAM3D::CellData::Phases),
m_GoodVoxelsArrayPath(DREAM3D::Defaults::ImageDataContainerName, DREAM3D::Defaults::CellAttributeMatrixName, DREAM3D::CellData::Mask),
m_CrystalStructuresArrayPath(DREAM3D::Defaults::ImageDataContainerName, DREAM3D::Defaults::CellEnsembleAttributeMatrixName, DREAM3D::EnsembleData::CrystalStructures),
m_Quats(NULL),
m_CellPhases(NULL),
m_GoodVoxels(NULL),
m_CrystalStructures(NULL)
{
	m_Seed = QDateTime::currentMSecsSinceEpoch();

	m_OrientationOps = SpaceGroupOps::getOrientationOpsQVector();

	featurecounts = NULL;

	// only setting up the child parameters because the parent constructor has already been called
	setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignmentMutualInformation::~AdaptiveAlignmentMutualInformation()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::setupFilterParameters()
{
	// getting the current parameters that were set by the parent and adding to it before resetting it
	FilterParameterVector parameters = getFilterParameters();
	parameters.push_front(DoubleFilterParameter::New("Misorientation Tolerance", "MisorientationTolerance", getMisorientationTolerance(), FilterParameter::Parameter));
	QStringList linkedProps("GoodVoxelsArrayPath");
	parameters.push_back(LinkedBooleanFilterParameter::New("Use Mask Array", "UseGoodVoxels", getUseGoodVoxels(), linkedProps, FilterParameter::Parameter));
	parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
	{
		DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(DREAM3D::TypeNames::Float, 4, DREAM3D::AttributeMatrixType::Cell, DREAM3D::GeometryType::ImageGeometry);
		parameters.push_back(DataArraySelectionFilterParameter::New("Quaternions", "QuatsArrayPath", getQuatsArrayPath(), FilterParameter::RequiredArray, req));
	}
  {
	  DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(DREAM3D::TypeNames::Int32, 1, DREAM3D::AttributeMatrixType::Cell, DREAM3D::GeometryType::ImageGeometry);
	  parameters.push_back(DataArraySelectionFilterParameter::New("Phases", "CellPhasesArrayPath", getCellPhasesArrayPath(), FilterParameter::RequiredArray, req));
  }
  {
	  DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(DREAM3D::TypeNames::Bool, 1, DREAM3D::AttributeMatrixType::Cell, DREAM3D::GeometryType::ImageGeometry);
	  parameters.push_back(DataArraySelectionFilterParameter::New("Mask", "GoodVoxelsArrayPath", getGoodVoxelsArrayPath(), FilterParameter::RequiredArray, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Cell Ensemble Data", FilterParameter::RequiredArray));
  {
	  DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(DREAM3D::TypeNames::UInt32, 1, DREAM3D::AttributeMatrixType::CellEnsemble, DREAM3D::GeometryType::ImageGeometry);

	  parameters.push_back(DataArraySelectionFilterParameter::New("Crystal Structures", "CrystalStructuresArrayPath", getCrystalStructuresArrayPath(), FilterParameter::RequiredArray, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
	AdaptiveAlignment::readFilterParameters(reader, index);
	reader->openFilterGroup(this, index);
	setCrystalStructuresArrayPath(reader->readDataArrayPath("CrystalStructuresArrayPath", getCrystalStructuresArrayPath()));
	setUseGoodVoxels(reader->readValue("UseGoodVoxels", getUseGoodVoxels()));
	setGoodVoxelsArrayPath(reader->readDataArrayPath("GoodVoxelsArrayPath", getGoodVoxelsArrayPath()));
	setCellPhasesArrayPath(reader->readDataArrayPath("CellPhasesArrayPath", getCellPhasesArrayPath()));
	setQuatsArrayPath(reader->readDataArrayPath("QuatsArrayPath", getQuatsArrayPath()));
	setMisorientationTolerance(reader->readValue("MisorientationTolerance", getMisorientationTolerance()));
	reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int AdaptiveAlignmentMutualInformation::writeFilterParameters(AbstractFilterParametersWriter* writer, int index)
{
	AdaptiveAlignment::writeFilterParameters(writer, index);
	writer->openFilterGroup(this, index);
	SIMPL_FILTER_WRITE_PARAMETER(CrystalStructuresArrayPath)
		SIMPL_FILTER_WRITE_PARAMETER(GoodVoxelsArrayPath)
		SIMPL_FILTER_WRITE_PARAMETER(UseGoodVoxels)
		SIMPL_FILTER_WRITE_PARAMETER(CellPhasesArrayPath)
		SIMPL_FILTER_WRITE_PARAMETER(QuatsArrayPath)
		SIMPL_FILTER_WRITE_PARAMETER(MisorientationTolerance)
		writer->closeFilterGroup();
	return ++index; // we want to return the next index that was just written to
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::dataCheck()
{
	setErrorCondition(0);

	// Set the DataContainerName and AttributematrixName for the Parent Class (AdaptiveAlignment) to Use.
	setDataContainerName(m_QuatsArrayPath.getDataContainerName());
	setCellAttributeMatrixName(m_QuatsArrayPath.getAttributeMatrixName());

	AdaptiveAlignment::dataCheck();
	if (getErrorCondition() < 0) { return; }

	INIT_DataArray(m_FeatureCounts, int32_t);

	QVector<DataArrayPath> dataArrayPaths;

	QVector<size_t> cDims(1, 4);
	m_QuatsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>, AbstractFilter>(this, getQuatsArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
	if (NULL != m_QuatsPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
	{
		m_Quats = m_QuatsPtr.lock()->getPointer(0);
	} /* Now assign the raw pointer to data from the DataArray<T> object */
	if (getErrorCondition() >= 0) { dataArrayPaths.push_back(getQuatsArrayPath()); }

	cDims[0] = 1;
	m_CellPhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray< int32_t>, AbstractFilter>(this, getCellPhasesArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
	if (NULL != m_CellPhasesPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
	{
		m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
	} /* Now assign the raw pointer to data from the DataArray<T> object */
	if (getErrorCondition() >= 0) { dataArrayPaths.push_back(getCellPhasesArrayPath()); }
	if (m_UseGoodVoxels == true)
	{
		m_GoodVoxelsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getGoodVoxelsArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
		if (NULL != m_GoodVoxelsPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
		{
			m_GoodVoxels = m_GoodVoxelsPtr.lock()->getPointer(0);
		} /* Now assign the raw pointer to data from the DataArray<T> object */
		if (getErrorCondition() >= 0) { dataArrayPaths.push_back(getGoodVoxelsArrayPath()); }
	}

	m_CrystalStructuresPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<unsigned int>, AbstractFilter>(this, getCrystalStructuresArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
	if (NULL != m_CrystalStructuresPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
	{
		m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
	} /* Now assign the raw pointer to data from the DataArray<T> object */

	getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::preflight()
{
	setInPreflight(true);
	emit preflightAboutToExecute();
	emit updateFilterParameters(this);
	dataCheck();
	emit preflightExecuted();
	setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts)
{
	DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());

	int64_t totalPoints = m->getAttributeMatrix(getCellAttributeMatrixName())->getNumTuples();
	m_MIFeaturesPtr = Int32ArrayType::CreateArray((totalPoints * 1), "_INTERNAL_USE_ONLY_MIFeatureIds");
	m_MIFeaturesPtr->initializeWithZeros();
	int32_t* miFeatureIds = m_MIFeaturesPtr->getPointer(0);

	size_t udims[3] = { 0, 0, 0 };
	m->getGeometryAs<ImageGeom>()->getDimensions(udims);

	uint64_t dims[3] =
	{
		static_cast<uint64_t>(udims[0]),
		static_cast<uint64_t>(udims[1]),
		static_cast<uint64_t>(udims[2]),
	};

	uint64_t maxstoredshifts = 1;
	if (xneedshifts.size() > 0) maxstoredshifts = 20;

	float disorientation = 0.0f;

	std::vector<std::vector<int64_t>>  newxshift(dims[2]);
	std::vector<std::vector<int64_t>>  newyshift(dims[2]);
	std::vector<std::vector<float>>  mindisorientation(dims[2]);
	for (uint64_t a = 1; a < dims[2]; a++)
	{
		newxshift[a].resize(maxstoredshifts, 0);
		newyshift[a].resize(maxstoredshifts, 0);
		mindisorientation[a].resize(maxstoredshifts, std::numeric_limits<float>::max());
	}

	float** mutualinfo12 = NULL;
	float* mutualinfo1 = NULL;
	float* mutualinfo2 = NULL;
	int32_t featurecount1 = 0, featurecount2 = 0;
	int64_t oldxshift = 0;
	int64_t oldyshift = 0;
	float count = 0.0f;
	uint64_t slice = 0;

	int32_t refgnum = 0, curgnum = 0;
	uint64_t refposition = 0;
	uint64_t curposition = 0;

	form_features_sections();

	// Allocate a 2D Array which will be reused from slice to slice
	// second dimension is assigned in each cycle separately
	std::vector<std::vector<bool> >  misorients(dims[0]);

	const uint64_t halfDim0 = static_cast<uint64_t>(dims[0] * 0.5f);
	const uint64_t halfDim1 = static_cast<uint64_t>(dims[1] * 0.5f);
	uint64_t progInt = 0;

	for (uint64_t iter = 1; iter < dims[2]; iter++)
	{
		progInt = static_cast<uint64_t>(iter * 100 / static_cast<float>(dims[2]));
		QString ss = QObject::tr("Aligning Anisotropic Sections || Determining Shifts || %1% Complete").arg(progInt);
		notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

		slice = (dims[2] - 1) - iter;
		featurecount1 = featurecounts[slice];
		featurecount2 = featurecounts[slice + 1];
		mutualinfo12 = new float *[featurecount1];
		mutualinfo1 = new float[featurecount1];
		mutualinfo2 = new float[featurecount2];

		for (int32_t a = 0; a < featurecount1; a++)
		{
			mutualinfo1[a] = 0.0f;
			mutualinfo12[a] = new float[featurecount2];
			for (int32_t b = 0; b < featurecount2; b++)
			{
				mutualinfo12[a][b] = 0.0f;
				mutualinfo2[b] = 0.0f;
			}
		}
		oldxshift = -1;
		oldyshift = -1;

		for (uint64_t i = 0; i < dims[0]; i++)
		{
			misorients[i].assign(dims[1], false);
		}

		while (newxshift[iter][0] != oldxshift || newyshift[iter][0] != oldyshift)
		{
			oldxshift = newxshift[iter][0];
			oldyshift = newyshift[iter][0];

			for (int32_t j = -3; j < 4; j++)
			{
				for (int32_t k = -3; k < 4; k++)
				{
					disorientation = 0;
					count = 0;
					if (llabs(k + oldxshift) < halfDim0 && llabs(j + oldyshift) < halfDim1 && misorients[k + oldxshift + halfDim0][j + oldyshift + halfDim1] == false)
					{
						for (uint64_t l = 0; l < dims[1]; l = l + 4)
						{
							for (uint64_t n = 0; n < dims[0]; n = n + 4)
							{
								if ((l + j + oldyshift) >= 0 && (l + j + oldyshift) < dims[1] && (n + k + oldxshift) >= 0 && (n + k + oldxshift) < dims[0])
								{
									refposition = ((slice + 1) * dims[0] * dims[1]) + (l * dims[0]) + n;
									curposition = (slice * dims[0] * dims[1]) + ((l + j + oldyshift) * dims[0]) + (n + k + oldxshift);
									refgnum = miFeatureIds[refposition];
									curgnum = miFeatureIds[curposition];
									if (curgnum >= 0 && refgnum >= 0)
									{
										mutualinfo12[curgnum][refgnum]++;
										mutualinfo1[curgnum]++;
										mutualinfo2[refgnum]++;
										count++;
									}
								}
								else
								{
									mutualinfo12[0][0]++;
									mutualinfo1[0]++;
									mutualinfo2[0]++;
								}
							}
						}
						float ha = 0.0f;
						float hb = 0.0f;
						float hab = 0.0f;
						for (int32_t b = 0; b < featurecount1; b++)
						{
							mutualinfo1[b] = mutualinfo1[b] / count;
							if (mutualinfo1[b] != 0) { ha = ha + mutualinfo1[b] * logf(mutualinfo1[b]); }
						}
						for (int32_t c = 0; c < featurecount2; c++)
						{
							mutualinfo2[c] = mutualinfo2[c] / float(count);
							if (mutualinfo2[c] != 0) { hb = hb + mutualinfo2[c] * logf(mutualinfo2[c]); }
						}
						for (int32_t b = 0; b < featurecount1; b++)
						{
							for (int32_t c = 0; c < featurecount2; c++)
							{
								mutualinfo12[b][c] = mutualinfo12[b][c] / count;
								if (mutualinfo12[b][c] != 0) { hab = hab + mutualinfo12[b][c] * logf(mutualinfo12[b][c]); }
								float value = 0.0f;
								if (mutualinfo1[b] > 0 && mutualinfo2[c] > 0) { value = (mutualinfo12[b][c] / (mutualinfo1[b] * mutualinfo2[c])); }
								if (value != 0) { disorientation = disorientation + (mutualinfo12[b][c] * logf(value)); }
							}
						}
						for (int32_t b = 0; b < featurecount1; b++)
						{
							for (int32_t c = 0; c < featurecount2; c++)
							{
								mutualinfo12[b][c] = 0.0f;
								mutualinfo1[b] = 0.0f;
								mutualinfo2[c] = 0.0f;
							}
						}
						disorientation = 1.0f / disorientation;
						misorients[k + oldxshift + halfDim0][j + oldyshift + halfDim1] = true;

						// compare the new shift with currently stored ones
						int64_t s = maxstoredshifts;
						while (s - 1 >= 0 && disorientation < mindisorientation[iter][s - 1])
						{
							s--;
						}

						// new shift is stored with index 's' in the arrays
						if (s < maxstoredshifts)
						{
							// lag the shifts already stored
							for (int64_t t = maxstoredshifts - 1; t > s; t--)
							{
								newxshift[iter][t] = newxshift[iter][t - 1];
								newyshift[iter][t] = newyshift[iter][t - 1];
								mindisorientation[iter][t] = mindisorientation[iter][t - 1];
							}
							// store the new shift
							newxshift[iter][s] = k + oldxshift;
							newyshift[iter][s] = j + oldyshift;
							mindisorientation[iter][s] = disorientation;
						}
					}
				}
			}
		}
		xshifts[iter] = xshifts[iter - 1] + newxshift[iter][0];
		yshifts[iter] = yshifts[iter - 1] + newyshift[iter][0];

		delete[] mutualinfo1;
		delete[] mutualinfo2;
		for (int32_t i = 0; i < featurecount1; i++)
		{
			delete mutualinfo12[i];
		}
		delete[] mutualinfo12;
		mutualinfo1 = NULL;
		mutualinfo2 = NULL;
		mutualinfo12 = NULL;
	}

	std::vector<uint64_t> curindex(dims[2], 0);

	// find corrected shifts
	if (xneedshifts.size() > 0)
	{
		QString ss = QObject::tr("Aligning Anisotropic Sections || Correcting shifts");
		notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

		std::vector<float> changedisorientation(dims[2], 0);
		std::vector<uint64_t> changeindex(dims[2], 0);
		std::vector<float> changeerror(dims[2], 0);

		std::vector<float> xshiftsest;	// cumulative x-shifts estimated from SEM images
		std::vector<float> yshiftsest;	// cumulative y-shifts estimated from SEM images

		float curerror = 0;
		float tolerance = 1.0f / static_cast<float>(dims[2]);

		// evaluate error between current shifts and desired shifts
		if (xneedshifts.size() == 1)      // error is computed as misagreement between slopes
		{
			curerror = compute_error1(dims[2], 0, xneedshifts[0], yneedshifts[0], newxshift, newyshift, curindex);
		}
		else if (xneedshifts.size() > 1)  // error is computed as misagreement with shifts estimated from SEM images
		{
			xshiftsest.resize(dims[2], 0);
			yshiftsest.resize(dims[2], 0);
			for (uint64_t iter = 1; iter < dims[2]; iter++)
			{
				xshiftsest[iter] = xshiftsest[iter - 1] + xneedshifts[iter - 1];
				yshiftsest[iter] = yshiftsest[iter - 1] + yneedshifts[iter - 1];
			}
			curerror = compute_error2(dims[2], 0, xshiftsest, yshiftsest, newxshift, newyshift, curindex);
		}

		// iterative selection of a candidate shift, recomputing of current candidates, evaluation of error
		if (curerror > tolerance)
		{
			float minchangedisorientation = 0;
			float minchangeerror = 0;
			int64_t minchangeindex = 0;
			int64_t minchangeiter = 0;
			float olderror = 0;
			float newerror = 0;
			uint64_t progInt = 0;

			do
			{
				QString ss = QObject::tr("Aligning Anisotropic Sections || Correcting shifts || Iteration %1").arg(++progInt);;
				notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
				if (getCancel() == true)
				{
					return;
				}

				olderror = curerror;
				for (uint64_t iter = 1; iter < dims[2]; iter++)
				{
					float newminerror = std::numeric_limits<float>::max();
					float newmindisorientation = std::numeric_limits<float>::max();
					uint64_t newminindex = 0;
					for (uint64_t index = curindex[iter] + 1; index < maxstoredshifts; index++)
					{
						// recompute error for the configuration with this candidate changed
						if (xneedshifts.size() == 1)
						{
							newerror = compute_error1(iter, index, xneedshifts[0], yneedshifts[0], newxshift, newyshift, curindex);
						}
						else if (xneedshifts.size() > 1)
						{
							newerror = compute_error2(iter, index, xshiftsest, yshiftsest, newxshift, newyshift, curindex);
						}

						// compare the new error with the best current error
						if (newerror < curerror &&
							mindisorientation[iter][index] / mindisorientation[iter][0] < newmindisorientation)
						{
							newminerror = newerror;
							newminindex = index;
							newmindisorientation = mindisorientation[iter][index] / mindisorientation[iter][0];
						}
					}
					// assign best error, corresponding index and disorientation value for this slice
					changeerror[iter] = newminerror;
					changeindex[iter] = newminindex;
					changedisorientation[iter] = newmindisorientation;
				}

				// among all slices, find the best candidate (with minimum disorientation change)
				minchangedisorientation = std::numeric_limits<float>::max() - 1;
				minchangeerror = std::numeric_limits<float>::max();
				minchangeindex = 0;
				minchangeiter = 0;
				for (uint64_t iter = 1; iter < dims[2]; iter++)
				{
					if (changeerror[iter] < curerror &&
						(changedisorientation[iter] < minchangedisorientation ||
						(changedisorientation[iter] == minchangedisorientation &&
						llabs(newxshift[iter][changeindex[iter]]) + llabs(newyshift[iter][changeindex[iter]]) < llabs(newxshift[iter][minchangeindex]) + llabs(newyshift[iter][minchangeindex]))))
					{
						minchangeiter = iter;
						minchangeindex = changeindex[iter];
						minchangedisorientation = changedisorientation[iter];
						minchangeerror = changeerror[iter];
					}
				}

				if (minchangeerror < curerror && minchangeerror >= tolerance)
				{
					// assign the best candidate
					changedisorientation[minchangeiter] = minchangedisorientation;
					curindex[minchangeiter] = minchangeindex;
					// reassign current error
					curerror = minchangeerror;
				}

			} while (minchangedisorientation < std::numeric_limits<float>::max() - 1 && curerror < olderror && curerror > tolerance);
		}
	}


	if (getWriteAlignmentShifts() == true)
	{
		std::ofstream outFile;
		outFile.open(getAlignmentShiftFileName().toLatin1().data());
		for (uint64_t iter = 1; iter < dims[2]; iter++)
		{
			slice = (dims[2] - 1) - iter;
			xshifts[iter] = xshifts[iter - 1] + newxshift[iter][curindex[iter]];
			yshifts[iter] = yshifts[iter - 1] + newyshift[iter][curindex[iter]];
			outFile << slice << "	" << slice + 1 << "	" << newxshift[iter][curindex[iter]] << "	" << newyshift[iter][curindex[iter]] << "	" << xshifts[iter] << "	" << yshifts[iter] << "\n";
		}
		outFile.close();
	}

	m->getAttributeMatrix(getCellAttributeMatrixName())->removeAttributeArray(DREAM3D::CellData::FeatureIds);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// find the error if the current shifts are changed by modification provided on slice 'iter' 
// where currently selected shift is substituted by the shift stored in position 'index'
float AdaptiveAlignmentMutualInformation::compute_error1(uint64_t iter, uint64_t index, float xneedtrend, float yneedtrend,
	std::vector<std::vector<int64_t>>& newxshift, std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex)
{
	int64_t n = curindex.size() - 1;

	int64_t xshifts = 0;
	int64_t yshifts = 0;

	double sumX = 0.0f;
	double sumX_2 = 0.0f;
	double x_sumY = 0.0f;
	double x_sumXY = 0.0f;
	double y_sumY = 0.0f;
	double y_sumXY = 0.0f;

	for (int64_t i = 1; i <= n; i++)
	{
		if (i != iter) // shifts without modification
		{
			xshifts += newxshift[i][curindex[i]];
			yshifts += newyshift[i][curindex[i]];
		}
		else           // modified shift
		{
			xshifts += newxshift[i][index];
			yshifts += newyshift[i][index];
		}
		// all values are divided by n for better numerical stability
		sumX += static_cast<double>(i) / static_cast<double>(n);
		sumX_2 += static_cast<double>(i * i) / static_cast<double>(n * n);
		x_sumY += static_cast<double>(xshifts) / static_cast<double>(n);
		x_sumXY += static_cast<double>(i * xshifts) / static_cast<double>(n * n);
		y_sumY += static_cast<double>(yshifts) / static_cast<double>(n);
		y_sumXY += static_cast<double>(i * yshifts) / static_cast<double>(n * n);
	}

	double xtrend = static_cast<double>((n * x_sumXY - sumX * x_sumY) / (n * sumX_2 - sumX * sumX));
	double ytrend = static_cast<double>((n * y_sumXY - sumX * y_sumY) / (n * sumX_2 - sumX * sumX));

	// error is computed from angular deviation as |tg(ang2 - ang1)| = |(tg(ang2) - tg(ang1)) / (1 + tg(ang1) * tg(ang2)|
	// where tg(ang1) = xneedtrend (resp. yneedtrend) and tg(ang2) = xtrend (resp. ytrend)
	return std::fabs((xtrend - xneedtrend) / (1 + xtrend * xneedtrend)) + std::fabs((ytrend - yneedtrend) / (1 + ytrend * yneedtrend));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float AdaptiveAlignmentMutualInformation::compute_error2(uint64_t iter, uint64_t index, std::vector<float>& xshiftsest, std::vector<float>& yshiftsest,
	std::vector<std::vector<int64_t>>& newxshift, std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex)
{
	uint64_t n = curindex.size() - 1;

	int64_t xshifts = 0;
	int64_t yshifts = 0;

	float error = 0;
	float xdif = 0;
	float ydif = 0;

	for (uint64_t i = 1; i <= n; i++)
	{
		if (i != iter)
		{
			xshifts += newxshift[i][curindex[i]];
			yshifts += newyshift[i][curindex[i]];
		}
		else
		{
			xshifts += newxshift[i][index];
			yshifts += newyshift[i][index];
		}

		xdif = static_cast<float>(xshifts)-xshiftsest[i];
		ydif = static_cast<float>(yshifts)-yshiftsest[i];
		error += (xdif * xdif + ydif * ydif);
	}
	return error;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::form_features_sections()
{
	SIMPL_RANDOMNG_NEW()
		DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());

	size_t udims[3] = { 0, 0, 0 };
	m->getGeometryAs<ImageGeom>()->getDimensions(udims);
	int64_t dims[3] =
	{
		static_cast<int64_t>(udims[0]),
		static_cast<int64_t>(udims[1]),
		static_cast<int64_t>(udims[2]),
	};

	int64_t point = 0;
	int64_t seed = 0;
	bool noseeds = false;
	int32_t featurecount = 1;
	int64_t neighbor = 0;
	QuatF q1 = QuaternionMathF::New();
	QuatF q2 = QuaternionMathF::New();
	QuatF* quats = reinterpret_cast<QuatF*>(m_Quats);
	float w = 0.0f;
	float n1 = 0.0f;
	float n2 = 0.0f;
	float n3 = 0.0f;
	int64_t randx = 0;
	int64_t randy = 0;
	bool good = false;
	int64_t x = 0, y = 0, z = 0;
	int64_t col = 0, row = 0;
	size_t size = 0;
	size_t initialVoxelsListSize = 1000;

	m_FeatureCounts->resize(dims[2]);
	featurecounts = m_FeatureCounts->getPointer(0);

	int32_t* miFeatureIds = m_MIFeaturesPtr->getPointer(0);

	std::vector<int64_t> voxelslist(initialVoxelsListSize, -1);
	int64_t neighpoints[4] = { 0, 0, 0, 0 };
	neighpoints[0] = -dims[0];
	neighpoints[1] = -1;
	neighpoints[2] = 1;
	neighpoints[3] = dims[0];

	uint32_t phase1 = 0, phase2 = 0;

	for (int64_t slice = 0; slice < dims[2]; slice++)
	{
		float prog = ((float)slice / dims[2]) * 100;
		QString ss = QObject::tr("Aligning Sections || Identifying Features on Sections || %1% Complete").arg(QString::number(prog, 'f', 0));
		notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
		featurecount = 1;
		noseeds = false;
		while (noseeds == false)
		{
			seed = -1;
			randx = int64_t(float(rg.genrand_res53()) * float(dims[0]));
			randy = int64_t(float(rg.genrand_res53()) * float(dims[1]));
			for (int64_t j = 0; j < dims[1]; ++j)
			{
				for (int64_t i = 0; i < dims[0]; ++i)
				{
					x = randx + i;
					y = randy + j;
					z = slice;
					if (x > dims[0] - 1) { x = x - dims[0]; }
					if (y > dims[1] - 1) { y = y - dims[1]; }
					point = (z * dims[0] * dims[1]) + (y * dims[0]) + x;
					if ((m_UseGoodVoxels == false || m_GoodVoxels[point] == true) && miFeatureIds[point] == 0 && m_CellPhases[point] > 0)
					{
						seed = point;
					}
					if (seed > -1) { break; }
				}
				if (seed > -1) { break; }
			}
			if (seed == -1) { noseeds = true; }
			if (seed >= 0)
			{
				size = 0;
				miFeatureIds[seed] = featurecount;
				voxelslist[size] = seed;
				size++;
				for (size_t j = 0; j < size; ++j)
				{
					int64_t currentpoint = voxelslist[j];
					col = currentpoint % dims[0];
					row = (currentpoint / dims[0]) % dims[1];
					QuaternionMathF::Copy(quats[currentpoint], q1);
					phase1 = m_CrystalStructures[m_CellPhases[currentpoint]];
					for (int32_t i = 0; i < 4; i++)
					{
						good = true;
						neighbor = currentpoint + neighpoints[i];
						if ((i == 0) && row == 0) { good = false; }
						if ((i == 3) && row == (dims[1] - 1)) { good = false; }
						if ((i == 1) && col == 0) { good = false; }
						if ((i == 2) && col == (dims[0] - 1)) { good = false; }
						if (good == true && miFeatureIds[neighbor] <= 0 && m_CellPhases[neighbor] > 0)
						{
							w = std::numeric_limits<float>::max();
							QuaternionMathF::Copy(quats[neighbor], q2);
							phase2 = m_CrystalStructures[m_CellPhases[neighbor]];
							if (phase1 == phase2)
							{
								w = m_OrientationOps[phase1]->getMisoQuat(q1, q2, n1, n2, n3);
							}
							if (w < m_MisorientationTolerance)
							{
								miFeatureIds[neighbor] = featurecount;
								voxelslist[size] = neighbor;
								size++;
								if (std::vector<int64_t>::size_type(size) >= voxelslist.size())
								{
									size = voxelslist.size();
									voxelslist.resize(size + initialVoxelsListSize);
									for (std::vector<int64_t>::size_type v = size; v < voxelslist.size(); ++v) { voxelslist[v] = -1; }
								}
							}
						}
					}
				}
				voxelslist.erase(std::remove(voxelslist.begin(), voxelslist.end(), -1), voxelslist.end());
				featurecount++;
				voxelslist.assign(initialVoxelsListSize, -1);
			}
		}
		featurecounts[slice] = featurecount;
	}
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMutualInformation::execute()
{
	setErrorCondition(0);
	dataCheck();
	if (getErrorCondition() < 0) { return; }

	// Converting the user defined tolerance to radians.
	m_MisorientationTolerance = m_MisorientationTolerance * SIMPLib::Constants::k_Pi / 180.0f;

	AdaptiveAlignment::execute();

	// If there is an error set this to something negative and also set a message
	notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer AdaptiveAlignmentMutualInformation::newFilterInstance(bool copyFilterParameters)
{
	AdaptiveAlignmentMutualInformation::Pointer filter = AdaptiveAlignmentMutualInformation::New();
	if (true == copyFilterParameters)
	{
		copyFilterParameterInstanceVariables(filter.get());
	}
	return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMutualInformation::getCompiledLibraryName()
{
	return AnisotropyConstants::AnisotropyBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMutualInformation::getBrandingString()
{
	return "Anisotropy";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMutualInformation::getFilterVersion()
{
	QString version;
	QTextStream vStream(&version);
	vStream << SIMPLib::Version::Major() << "." << SIMPLib::Version::Minor() << "." << SIMPLib::Version::Patch();
	return version;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMutualInformation::getGroupName()
{
	return DREAM3D::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMutualInformation::getSubGroupName()
{
	return "Anisotropy";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMutualInformation::getHumanLabel()
{
	return "Adaptive Alignment (Mutual Information)";
}
