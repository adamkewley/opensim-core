#ifndef OPENSIM_SIMULATION_UTILITIES_H_
#define OPENSIM_SIMULATION_UTILITIES_H_
/* -------------------------------------------------------------------------- *
 *                     OpenSim:  SimulationUtilities.h                        *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2018 Stanford University and the Authors                *
 * Author(s): OpenSim Team                                                    *
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

#include "Model/Model.h"
#include "Manager/Manager.h"
#include <simbody/internal/Visualizer_InputListener.h>

namespace OpenSim {

/// @name General-purpose simulation driver for OpenSim models
/// @{
/** Simulate a model from an initial state and return the final state.
    If the model's useVisualizer flag is true, the user is repeatedly prompted
    to either begin simulating or quit. The provided state is not updated but
    the final state is returned at the end of the simulation, when finalTime is
    reached. %Set saveStatesFile=true to save the states to a storage file as:
    "<model_name>_states.sto". */
inline SimTK::State simulate(Model& model,
    const SimTK::State& initialState,
    double finalTime,
    bool saveStatesFile = false)
{
    // Returned state begins as a copy of the initial state
    SimTK::State state = initialState;
    SimTK::Visualizer::InputSilo* silo;

    bool simulateOnce = true;

    // Ensure the final time is in the future.
    const double initialTime = initialState.getTime();
    if (finalTime <= initialTime) {
        std::cout << "The final time must be in the future (current time is "
                  << initialTime << "); simulation aborted." << std::endl;
        return state;
    }

    // Configure the visualizer.
    if (model.getUseVisualizer()) {
        SimTK::Visualizer& viz = model.updVisualizer().updSimbodyVisualizer();
        // We use the input silo to get key presses.
        silo = &model.updVisualizer().updInputSilo();

        SimTK::DecorativeText help("Press any key to start a new simulation; "
            "ESC to quit.");
        help.setIsScreenText(true);
        viz.addDecoration(SimTK::MobilizedBodyIndex(0), SimTK::Vec3(0), help);

        viz.setShowSimTime(true);
        viz.drawFrameNow(state);
        std::cout << "A visualizer window has opened." << std::endl;

        // if visualizing enable replay
        simulateOnce = false;
    }

    // Simulate until the user presses ESC (or enters 'q' if visualization has
    // been disabled).
    do {
        if (model.getUseVisualizer()) {
            // Get a key press.
            silo->clear(); // Ignore any previous key presses.
            unsigned key, modifiers;
            silo->waitForKeyHit(key, modifiers);
            if (key == SimTK::Visualizer::InputListener::KeyEsc) { break; }
        }

        // reset the state to the initial state
        state = initialState;
        // Set up manager and simulate.
        SimTK::RungeKuttaMersonIntegrator integrator(model.getSystem());
        Manager manager(model, integrator);
        state.setTime(initialTime);
        manager.initialize(state);
        state = manager.integrate(finalTime);

        // Save the states to a storage file (if requested).
        if (saveStatesFile) {
            manager.getStateStorage().print(model.getName() + "_states.sto");
        }
    } while (!simulateOnce);

    return state;
}
/// @}

/** This function can be used to upgrade MOT files generated with versions
    3.3 and earlier in which some data columns are associated with coordinates
    that were incorrectly marked as Rotational (rather than Coupled). Specific
    instances of the issue are the patella coordinate in the Rajagopal 2015 and
    leg6dof9musc models. In these cases, the patella will visualize incorrectly
    in the GUI when replaying the kinematics from the MOT file, and static 
    optimization results.
    Moving from OpenSim 3.3 (and earlier) to OpenSim 4.0, MotionTypes for
    Coordinates are now determined strictly by the coordinates' owning Joint.
    In older models, the MotionType, particularly for CustomJoints, were user-
    specified. That entailed in some cases, incorrectly labeling a Coordinate
    as being Rotational, for example, when it is in fact Coupled. For the above
    models, for example, the patella Coordinate had been user-specified to be
    Rotational, but the angle of the patella about the Z-axis of the patella
    body, is a spline function (e.g. coupled function) of the knee_patella
    Coordinate. Thus, knee_patella Coordinate does not represent a Cartesian
    rotational coordinate and is not classified as Rotational. Use this utility
    to remove any unit conversions from Coordinates that were incorrectly labeled
    as Rotational in the past. For these Coordinates only, the utility will undo
    the incorrect radians to degrees conversion. */
inline void updateKinematicsFilesForUpdatedModel(const Model& model, 
                const std::vector<std::string>& filePaths, std::string suffix="")
{
    // Cycle through the data files 
    for (auto filePath : filePaths) {
        Storage motion(filePath);
        Storage* updatedMotion = 
            updateKinematicsStorageForUpdatedModel(model, motion);

        if (updatedMotion == nullptr) {
            continue; //no update was required, move on to next file
        }

        std::string outFilePath = filePath;
        if (suffix.size()) {
            auto back = filePath.rfind(".");
            outFilePath = filePath.substr(0, back - 1) + suffix + 
                            filePath.substr(back+1);
        }
        std::cout << "Writing converted motion '" << filePath << "' to '"
            << outFilePath << "'." << std::endl;

        updatedMotion->print(outFilePath);
    }
}

inline Storage* updateKinematicsStorageForUpdatedModel(const Model& model,
                                                      const Storage &kinematics)
{
    // There is no issue if the kinematics are in internal values (i.e. not 
    // converted to degrees
    if(!kinematics.isInDegrees())
    if (model.getDocumentFileVersion() >= 30415) {
        throw Exception("updateKinematicsStorageForUpdatedModel has no updates "
            "to make because the model '" + model.getName() + "'is up-to-date.\n"
            "If input motion files were generated with this model version, there is "
            "nothing further to be done. Otherwise, provide the original model "
            "file used to generate the motion files and try again.");
    }

    std::vector<const Coordinate*> problemCoords;
    auto coordinates = model.getComponentList<Coordinate>();
    for (auto& coord : coordinates) {
        const Coordinate::MotionType oldMotionType =
            coord.getUserSpecifiedMotionTypePriorTo40();
        const Coordinate::MotionType motionType = coord.getMotionType();

        if ((oldMotionType != Coordinate::MotionType::Undefined) &&
            (oldMotionType != motionType)) {
            problemCoords.push_back(&coord);
        }
    }

    if (problemCoords.size() == 0)
        return nullptr;

    Storage* updatedKinematics = kinematics.clone();
    // Cycle the inconsistent Coordinates
    for (auto coord : problemCoords) {
        // Get the corresponding column of data and if in degrees
        // undo the radians to degrees conversion on that column.
        int ix = updatedKinematics->getStateIndex(coord->getName());

        if (ix < 0) {
            std::cout << "updateKinematicsStorageForUpdatedModel(): motion '"
                << kinematics.getName() << "' does not contain inconsistent "
                << "coordinate '" << coord->getName() << "'." << std::endl;
        } 
        else {
            // convert this column back to internal values by undoing the
            // 180/pi conversion to degrees
            updatedKinematics->multiplyColumn(ix, SimTK_DTR);
        }
    }

}

} // end of namespace OpenSim

#endif // OPENSIM_SIMULATION_UTILITIES_H_
