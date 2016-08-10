/* ============================================================================
* Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
* Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
*
* Redistribution and use in source and binary forms, with or without
* modification,
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

#ifndef _AdaptiveAlignmentFeature_H_
#define _AdaptiveAlignmentFeature_H_

#include "SIMPLib/Common/AbstractFilter.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/SIMPLib.h"

#include "Anisotropy/AnisotropyConstants.h"

#include "Anisotropy/AnisotropyFilters/AdaptiveAlignment.h"

/**
* @brief The AdaptiveAlignmentFeature class. See [Filter documentation](@ref AdaptiveAlignmentfeature) for details.
*/
class AdaptiveAlignmentFeature : public AdaptiveAlignment {
    Q_OBJECT
  public : SIMPL_SHARED_POINTERS(AdaptiveAlignmentFeature)
    SIMPL_STATIC_NEW_MACRO(AdaptiveAlignmentFeature)
    SIMPL_TYPE_MACRO_SUPER(AdaptiveAlignmentFeature,
                           AdaptiveAlignment)

    virtual ~AdaptiveAlignmentFeature();

    SIMPL_FILTER_PARAMETER(DataArrayPath, GoodVoxelsArrayPath)
    Q_PROPERTY(DataArrayPath GoodVoxelsArrayPath READ getGoodVoxelsArrayPath WRITE setGoodVoxelsArrayPath)

    /**
    * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
    */
    virtual const QString getCompiledLibraryName();

    /**
    * @brief getBrandingString Returns the branding string for the filter, which is
    * a tag
    * used to denote the filter's association with specific plugins
    * @return Branding string
    */
    virtual const QString getBrandingString();

    /**
    * @brief getFilterVersion Returns a version string for this filter. Default
    * value is an empty string.
    * @return
    */
    virtual const QString getFilterVersion();

    /**
    * @brief newFilterInstance Reimplemented from @see AbstractFilter class
    */
    virtual AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters);

    /**
    * @brief getGroupName Reimplemented from @see AbstractFilter class
    */
    virtual const QString getGroupName();

    /**
    * @brief getSubGroupName Reimplemented from @see AbstractFilter class
    */
    virtual const QString getSubGroupName();

    /**
    * @brief getHumanLabel Reimplemented from @see AbstractFilter class
    */
    virtual const QString getHumanLabel();

    /**
    * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
    */
    virtual void setupFilterParameters();

    /**
    * @brief readFilterParameters Reimplemented from @see AbstractFilter class
    */
    virtual void readFilterParameters(AbstractFilterParametersReader *reader, int index);

    /**
    * @brief execute Reimplemented from @see AbstractFilter class
    */
    virtual void execute();

    /**
    * @brief preflight Reimplemented from @see AbstractFilter class
    */
    virtual void preflight();

  protected:
    AdaptiveAlignmentFeature();

    /**
    * @brief dataCheck Checks for the appropriate parameter values and
    * availability of arrays
    */
    void dataCheck();

    /**
    * @brief Initializes all the private instance variables.
    */
    void initialize();

    /**
    * @brief find_shifts Reimplemented from @see AdaptiveAlignment class
    * @param xshifts
    * @param yshifts
    * @param xneedshifts
    * @param yneedshifts
    */
    virtual void find_shifts(std::vector<int64_t> &xshifts,
                             std::vector<int64_t> &yshifts,
                             std::vector<float> &xneedshifts,
                             std::vector<float> &yneedshifts);

    /**
    * @brief compute_error1 Determines error between the current and desired shifts
    * (based on discrepancy of slopes)
    * @param iter
    * @param index
    * @param xneedtrend
    * @param yneedtrend
    * @param newxshift
    * @param newyshift
    * @param curindex
    * @return
    */
    virtual float compute_error1(uint64_t iter, uint64_t index, float xneedtrend,
                                 float yneedtrend,
                                 std::vector<std::vector<int64_t>> &newxshift,
                                 std::vector<std::vector<int64_t>> &newyshift,
                                 std::vector<uint64_t> &curindex);

    /**
    * @brief compute_error1 Determines error between the current and desired shifts
    * (based on discrepancy from the estimates obtained from SEM images)
    * @param iter
    * @param index
    * @param xshiftsest
    * @param yshiftsest
    * @param newxshift
    * @param newyshift
    * @param curindex
    * @return
    */
    virtual float compute_error2(uint64_t iter, uint64_t index,
                                 std::vector<float> &xshiftsest,
                                 std::vector<float> &yshiftsest,
                                 std::vector<std::vector<int64_t>> &newxshift,
                                 std::vector<std::vector<int64_t>> &newyshift,
                                 std::vector<uint64_t> &curindex);

  private:
    DEFINE_DATAARRAY_VARIABLE(bool, GoodVoxels)

    AdaptiveAlignmentFeature(const AdaptiveAlignmentFeature &); // Copy Constructor Not Implemented
    void operator=(const AdaptiveAlignmentFeature &); // Operator '=' Not Implemented
};

#endif /* AdaptiveAlignmentFeature_H_ */
