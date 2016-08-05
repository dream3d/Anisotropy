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
* Neither the name of the Czech Academy of Sciences, nor the names of its
* contributors may be used to endorse or promote products derived from this
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

#include "AdaptiveAlignmentMisorientation.h"

#include <fstream>

#include <QtCore/QDateTime>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/SIMPLibVersion.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "Anisotropy/AnisotropyConstants.h"

#include "moc_AdaptiveAlignmentMisorientation.cpp"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignmentMisorientation::AdaptiveAlignmentMisorientation() :
AdaptiveAlignment(),
m_MisorientationTolerance(5.0f),
m_UseGoodVoxels(true),
m_QuatsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Quats),
m_CellPhasesArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Phases),
m_GoodVoxelsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Mask),
m_CrystalStructuresArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellEnsembleAttributeMatrixName, SIMPL::EnsembleData::CrystalStructures),
m_Quats(NULL),
m_CellPhases(NULL),
m_GoodVoxels(NULL),
m_CrystalStructures(NULL)
{
  m_RandomSeed = QDateTime::currentMSecsSinceEpoch();

  m_OrientationOps = SpaceGroupOps::getOrientationOpsQVector();

  // only setting up the child parameters because the parent constructor has already been called
  setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignmentMisorientation::~AdaptiveAlignmentMisorientation()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMisorientation::setupFilterParameters()
{
  // getting the current parameters that were set by the parent and adding to it before resetting it
  FilterParameterVector parameters = getFilterParameters();
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Misorientation Tolerance (Degrees)", MisorientationTolerance, FilterParameter::Parameter, AdaptiveAlignmentMisorientation));
  QStringList linkedProps("GoodVoxelsArrayPath");
  parameters.push_back(LinkedBooleanFilterParameter::New("Use Mask Array", "UseGoodVoxels", getUseGoodVoxels(), FilterParameter::Parameter, SIMPL_BIND_SETTER(AdaptiveAlignmentMisorientation, this, UseGoodVoxels), SIMPL_BIND_GETTER(AdaptiveAlignmentMisorientation, this, UseGoodVoxels), linkedProps));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 4, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Quaternions", QuatsArrayPath, FilterParameter::RequiredArray, AdaptiveAlignmentMisorientation, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", CellPhasesArrayPath, FilterParameter::RequiredArray, AdaptiveAlignmentMisorientation, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", GoodVoxelsArrayPath, FilterParameter::RequiredArray, AdaptiveAlignmentMisorientation, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Cell Ensemble Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::UInt32, 1, SIMPL::AttributeMatrixType::CellEnsemble, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Crystal Structures", CrystalStructuresArrayPath, FilterParameter::RequiredArray, AdaptiveAlignmentMisorientation, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMisorientation::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  AdaptiveAlignment::readFilterParameters(reader, index);
  reader->openFilterGroup(this, index);
  setCrystalStructuresArrayPath(reader->readDataArrayPath("CrystalStructuresArrayPath", getCrystalStructuresArrayPath()));
  setGoodVoxelsArrayPath(reader->readDataArrayPath("GoodVoxelsArrayPath", getGoodVoxelsArrayPath()));
  setUseGoodVoxels(reader->readValue("UseGoodVoxels", getUseGoodVoxels()));
  setCellPhasesArrayPath(reader->readDataArrayPath("CellPhasesArrayPath", getCellPhasesArrayPath()));
  setQuatsArrayPath(reader->readDataArrayPath("QuatsArrayPath", getQuatsArrayPath()));
  setMisorientationTolerance(reader->readValue("MisorientationTolerance", getMisorientationTolerance()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMisorientation::initialize()
{
  m_RandomSeed = QDateTime::currentMSecsSinceEpoch();

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMisorientation::dataCheck()
{
  setErrorCondition(0);
  initialize();

  // Set the DataContainerName and AttributematrixName for the Parent Class (AlignSections) to Use.
  setDataContainerName(m_QuatsArrayPath.getDataContainerName());
  setCellAttributeMatrixName(m_QuatsArrayPath.getAttributeMatrixName());

  AdaptiveAlignment::dataCheck();
  if (getErrorCondition() < 0) { return; }

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

  m_CrystalStructuresPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<uint32_t>, AbstractFilter>(this, getCrystalStructuresArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if (NULL != m_CrystalStructuresPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  {
    m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMisorientation::preflight()
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
void AdaptiveAlignmentMisorientation::find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());

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

  int64_t oldxshift = 0;
  int64_t oldyshift = 0;
  float count = 0.0f;
  uint64_t slice = 0;
  float w = 0.0f;
  float n1 = 0.0f, n2 = 0.0f, n3 = 0.0f;
  QuatF q1 = QuaternionMathF::New();
  QuatF q2 = QuaternionMathF::New();
  uint64_t refposition = 0;
  uint64_t curposition = 0;
  QuatF* quats = reinterpret_cast<QuatF*>(m_Quats);

  uint32_t phase1 = 0, phase2 = 0;
  uint64_t progInt = 0;

  // Allocate a 2D Array which will be reused from slice to slice
  // second dimension is assigned in each cycle separately
  std::vector<std::vector<bool> >  misorients(dims[0]);

  const uint64_t halfDim0 = static_cast<uint64_t>(dims[0] * 0.5f);
  const uint64_t halfDim1 = static_cast<uint64_t>(dims[1] * 0.5f);

  // Loop over the Z Direction
  for (uint64_t iter = 1; iter < dims[2]; iter++)
  {
    progInt = static_cast<uint64_t>(iter * 100 / static_cast<float>(dims[2]));
    QString ss = QObject::tr("Aligning Anisotropic Sections || Determining Shifts || %1% Complete").arg(progInt);
    notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
    if (getCancel() == true)
    {
      return;
    }

    slice = (dims[2] - 1) - iter;
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

      for (int32_t j = -3; j <= 3; j++)
      {
        for (int32_t k = -3; k <= 3; k++)
        {
          disorientation = 0.0f;
          count = 0.0f;
          if (llabs(k + oldxshift) < halfDim0 && llabs(j + oldyshift) < halfDim1 && misorients[k + oldxshift + halfDim0][j + oldyshift + halfDim1] == false)
          {
            for (uint64_t l = 0; l < dims[1]; l = l + 4)
            {
              for (uint64_t n = 0; n < dims[0]; n = n + 4)
              {
                if (int64_t((l + j + oldyshift)) >= 0 && (l + j + oldyshift) < dims[1] && int64_t((n + k + oldxshift)) >= 0 && (n + k + oldxshift) < dims[0])
                {
                  count++;
                  refposition = ((slice + 1) * dims[0] * dims[1]) + (l * dims[0]) + n;
                  curposition = (slice * dims[0] * dims[1]) + ((l + j + oldyshift) * dims[0]) + (n + k + oldxshift);
                  if (m_UseGoodVoxels == false || (m_GoodVoxels[refposition] == true && m_GoodVoxels[curposition] == true))
                  {
                    w = std::numeric_limits<float>::max();
                    if (m_CellPhases[refposition] > 0 && m_CellPhases[curposition] > 0)
                    {
                      QuaternionMathF::Copy(quats[refposition], q1);
                      phase1 = m_CrystalStructures[m_CellPhases[refposition]];
                      QuaternionMathF::Copy(quats[curposition], q2);
                      phase2 = m_CrystalStructures[m_CellPhases[curposition]];
                      if (phase1 == phase2 && phase1 < static_cast<uint32_t>(m_OrientationOps.size()))
                      {
                        w = m_OrientationOps[phase1]->getMisoQuat(q1, q2, n1, n2, n3);
                      }
                    }
                    if (w > m_MisorientationTolerance) { disorientation++; }
                  }
                  if (m_UseGoodVoxels == true)
                  {
                    if (m_GoodVoxels[refposition] == true && m_GoodVoxels[curposition] == false) { disorientation++; }
                    if (m_GoodVoxels[refposition] == false && m_GoodVoxels[curposition] == true) { disorientation++; }
                  }
                }
                else
                {

                }
              }
            }

            disorientation = disorientation / count;
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

    float curerror = 0.0f;
    float tolerance = 0.0f;

    // evaluate error between current shifts and desired shifts
    if (xneedshifts.size() == 1)      // error is computed as misagreement between slopes
    {
      tolerance = 1.0f / static_cast<float>(dims[2] - 1);
      curerror = compute_error1(dims[2], 0, xneedshifts[0], yneedshifts[0], newxshift, newyshift, curindex);
    }
    else if (xneedshifts.size() > 1)  // error is computed as misagreement with shifts estimated from SEM images
    {
      tolerance = 1.0f;
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
        QString ss = QObject::tr("Aligning Anisotropic Sections || Correcting Shifts || Iteration %1").arg(++progInt);;
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
              mindisorientation[iter][index] - mindisorientation[iter][0] < newmindisorientation)
            {
              newminerror = newerror;
              newminindex = index;
              newmindisorientation = mindisorientation[iter][index] - mindisorientation[iter][0];
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
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// find the error if the current shifts are changed by modification provided on slice 'iter'
// where currently selected shift is substituted by the shift stored in position 'index'
float AdaptiveAlignmentMisorientation::compute_error1(uint64_t iter, uint64_t index, float xneedtrend, float yneedtrend,
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

  // error is computed from angular deviation as
  // |xang1 - xang2| + |yang1 - yang2| = |arctg(xtrend) - arctg(xneedtrend)| + |arctg(ytrend) - arctg(yneedtrend)|
  return std::fabs(atan(xtrend) - atan(xneedtrend)) + std::fabs(atan(ytrend) - atan(yneedtrend));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float AdaptiveAlignmentMisorientation::compute_error2(uint64_t iter, uint64_t index, std::vector<float>& xshiftsest, std::vector<float>& yshiftsest,
  std::vector<std::vector<int64_t>>& newxshift, std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex)
{
  uint64_t n = curindex.size() - 1;

  int64_t xshifts = 0;
  int64_t yshifts = 0;

  float error = 0;
  float xdif = 0;
  float ydif = 0;
  float divide = static_cast<float> (2 * n);

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
    error += (xdif * xdif + ydif * ydif) / divide;
  }
  return error;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentMisorientation::execute()
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
AbstractFilter::Pointer AdaptiveAlignmentMisorientation::newFilterInstance(bool copyFilterParameters)
{
  AdaptiveAlignmentMisorientation::Pointer filter = AdaptiveAlignmentMisorientation::New();
  if (true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMisorientation::getCompiledLibraryName()
{
  return AnisotropyConstants::AnisotropyBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMisorientation::getBrandingString()
{
  return "Anisotropy";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMisorientation::getFilterVersion()
{
  QString version;
  QTextStream vStream(&version);
  vStream << SIMPLib::Version::Major() << "." << SIMPLib::Version::Minor() << "." << SIMPLib::Version::Patch();
  return version;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMisorientation::getGroupName()
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMisorientation::getSubGroupName()
{
  return AnisotropyConstants::FilterSubGroups::AnisotropicAlignment;;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentMisorientation::getHumanLabel()
{
  return "Adaptive Alignment (Misorientation)";
}
