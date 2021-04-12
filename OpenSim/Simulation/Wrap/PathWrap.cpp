/* -------------------------------------------------------------------------- *
 *                           OpenSim:  PathWrap.cpp                           *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan                                                      *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

//=============================================================================
// INCLUDES
//=============================================================================
#include "PathWrap.h"
#include <OpenSim/Simulation/Model/Model.h>

#include <array>

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;

static constexpr std::array<const char*, 3> wrapMethod2WrapName = { "hybrid", "midpoint", "axial" };
static_assert(static_cast<size_t>(OpenSim::PathWrap::hybrid) == 0);
static_assert(static_cast<size_t>(OpenSim::PathWrap::midpoint) == 1);
static_assert(static_cast<size_t>(OpenSim::PathWrap::axial) == 2);

static const char* wrapMethodToString(OpenSim::PathWrap::WrapMethod method) {
    SimTK_ASSERT(static_cast<size_t>(method) < wrapMethod2WrapName.size());
    return wrapMethod2WrapName[static_cast<size_t>(method)];
}

static void appendAcceptedWrapNames(std::ostringstream& ss) {
    const char* prefix = "'";
    for (const auto& name : wrapmethod2wrapname_strict) {
        ss << name << '\'';
        prefix = ", '";
    }
}

static const std::unordered_map<std::string, OpenSim::PathWrap::WrapMethod>
    wrapName2WrapMethod = {
        {"hybrid", OpenSim::PathWrap::hybrid},
        {"Hybrid", OpenSim::PathWrap::hybrid},
        {"HYBRID", OpenSim::PathWrap::hybrid},
        {"midpoint", OpenSim::PathWrap::midpoint},
        {"Midpoint", OpenSim::PathWrap::midpoint},
        {"MIDPOINT", OpenSim::PathWrap::midpoint},
        {"axial", OpenSim::PathWrap::axial},
        {"Axial", OpenSim::PathWrap::axial},
        {"AXIAL", OpenSim::PathWrap::axial},
    };

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor.
 */
PathWrap::PathWrap() : ModelComponent()
{
    setNull();
    constructProperties();
}

//_____________________________________________________________________________
/**
 * Destructor.
 */
PathWrap::~PathWrap()
{
}

//=============================================================================
// CONSTRUCTION METHODS
//=============================================================================
//_____________________________________________________________________________
/**
 * Set the data members of this PathWrap to their null values.
 */
void PathWrap::setNull()
{
    resetPreviousWrap();
}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void PathWrap::constructProperties()
{
    constructProperty_wrap_object("");
    constructProperty_method(wrapMethodToString(WrapMethod::hybrid));
    OpenSim::Array<int> range(-1, 2);
    constructProperty_range(range);
}

void PathWrap::extendConnectToModel(Model& model)
{
    Super::extendConnectToModel(model);

    _path = dynamic_cast<const GeometryPath*>(&getOwner());

    if (_path == nullptr) {
        std::stringstream ss;
        ss << "PathWrap '" << getAbsolutePathString()
           << "' must have a GeometryPath as its owner.";
        OPENSIM_THROW(Exception, std::move(ss).str());
    }

    // assign the wrap object from the specified wrap object name
    {
        std::string const& woName = getWrapObjectName();

        bool found = false;
        for (auto const& wo : model.getComponentList<OpenSim::WrapObject>()) {
            if (wo.getName() == woName) {
                _wrapObject = &wo;
                updWrapPoint1().setParentFrame(wo.getFrame());
                updWrapPoint1().setWrapObject(&wo);
                updWrapPoint2().setParentFrame(wo.getFrame());
                updWrapPoint2().setWrapObject(&wo);
                found = true;
                break;
            }
        }

        if (!found) {
            std::stringstream ss;
            ss << "Cannot connect PathWrap '" << getAbsolutePathString()
               << "' to the wrap object '" << woName
               << "': the specified wrap object could not be found in the model";
            OPENSIM_THROW(Exception, std::move(ss).str());
        }
    }

    // assign _method from the specified method name string
    {
        // edge-case: the method name isn't specified
        if (get_method().empty()) {
            std::stringstream ss;
            ss << "No method name specified for PathWrap '" << getAbsolutePathString()
               << "': must be one of: ";
            appendAcceptedWrapNames(ss);

            OPENSIM_THROW(Exception, std::move(ss).str());
        }

        // edge-case: "Unassigned" is specified: the implementation
        // should choose a reasonable default
        if (get_method() == "Unassigned") {
            upd_method() = wrapMethodToString(WrapMethod::hybrid);
        }

        // normal-case: something was specified: it should be in the LUT
        auto it = wrapName2WrapMethod.find(get_method());

        if (it == wrapName2WrapMethod.end()) {
            std::stringstream ss;
            ss << "The method name '" << get_method() << "' for PathWrap '"
               << getAbsolutePathString() << "' is invalid. Allowed values are: ";
            appendAcceptedWrapNames(ss);

            OPENSIM_THROW(Exception, std::move(ss).str());
        }

        _method = it->second;
    }
}

void PathWrap::setStartPoint( const SimTK::State& s, int aIndex)
{
    if ((aIndex != get_range(0)) && 
        (aIndex == -1 || get_range(1) == -1 || (aIndex >= 1 && aIndex <= get_range(1))))
    {
        upd_range(0) = aIndex;
    }
}

void PathWrap::setEndPoint( const SimTK::State& s, int aIndex)
{
    if ((aIndex != get_range(1)) && 
        (aIndex == -1 || get_range(0) == -1 || (aIndex >= get_range(0) && aIndex <= _path->getPathPointSet().getSize())))
    {
        upd_range(1) = aIndex;
    }
}

void PathWrap::resetPreviousWrap()
{
    _previousWrap.startPoint = -1;
    _previousWrap.endPoint = -1;

    _previousWrap.wrap_pts.setSize(0);
    _previousWrap.wrap_path_length = 0.0;

    int i;
    for (i = 0; i < 3; i++) {
        _previousWrap.r1[i] = -std::numeric_limits<SimTK::Real>::infinity();
        _previousWrap.r2[i] = -std::numeric_limits<SimTK::Real>::infinity();
        _previousWrap.sv[i] = -std::numeric_limits<SimTK::Real>::infinity();
    }
}

void PathWrap::setPreviousWrap(const WrapResult& aWrapResult)
{
    _previousWrap = aWrapResult;
}

void PathWrap::setWrapObject(WrapObject& aWrapObject)
{
    _wrapObject = &aWrapObject;
    upd_wrap_object() = aWrapObject.getName();
}

void PathWrap::setMethod(WrapMethod aMethod)
{
    _method = aMethod;
    upd_method() = wrapMethodToString(aMethod);
}
