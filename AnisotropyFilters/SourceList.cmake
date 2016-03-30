#--////////////////////////////////////////////////////////////////////////////
#-- Your License or copyright can go here
#--////////////////////////////////////////////////////////////////////////////

set(_filterGroupName AnisotropyFilters)
set(${_filterGroupName}_FILTERS_HDRS "")

START_FILTER_GROUP(${Anisotropy_BINARY_DIR} "${_filterGroupName}" "Anisotropy Filters")

#---------
# List your public filters here

set(_PublicFilters
  AdaptiveAlignmentFeature
  AdaptiveAlignmentMisorientation
  AdaptiveAlignmentMutualInformation
  SteinerCompact
)

#--------------
# Loop on all the filters adding each one. In this loop we default to making each filter exposed in the user
# interface in DREAM3D. If you want to have the filter compiled but NOT exposed to the user then use the next loop
foreach(f ${_PublicFilters} )
  ADD_DREAM3D_FILTER(  "Anisotropy" "Anisotropy"
                        ${_filterGroupName} ${f}
                        ${Anisotropy_SOURCE_DIR}/Documentation/${_filterGroupName}/${f}.md TRUE)
endforeach()


#---------------
# This is the list of Private Filters. These filters are available from other filters but the user will not
# be able to use them from the DREAM3D user interface.
set(_PrivateFilters
  AdaptiveAlignment
)

#-----------------
# Loop on the Private Filters adding each one to the DREAM3DLib project so that it gets compiled.
foreach(f ${_PrivateFilters} )
  ADD_DREAM3D_FILTER(  "Anisotropy" "Anisotropy"
                        ${_filterGroupName} ${f}
                        ${Anisotropy_SOURCE_DIR}/Documentation/${_filterGroupName}/${f}.md FALSE)
endforeach()

END_FILTER_GROUP(${Anisotropy_BINARY_DIR} "${_filterGroupName}" "Anisotropy Filters")

