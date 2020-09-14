#ifndef _WrapMath_h_
#define _WrapMath_h_
/* -------------------------------------------------------------------------- *
 *                            OpenSim:  WrapMath.h                            *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Frank C. Anderson                                               *
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

/* Note: This code was originally developed by Realistic Dynamics Inc. 
 * Author: Frank C. Anderson 
 */

#include <OpenSim/Simulation/osimSimulationDLL.h>
#include <SimTKcommon/SmallMatrix.h>


namespace OpenSim { 

/** @cond **/ // hide from Doxygen

//=============================================================================
//=============================================================================
/**
 * This class provides basic math functions and constants for wrapping surfaces.
 */
struct OSIMSIMULATION_API WrapMath {

    /**
     * Compute the intersection between a line (p1->p2) and another line
     * (p3->p4). If the lines do not intersect, this function returns the
     * closest point on each line to the other line.
     *
     * @param p1 first point on first line
     * @param p2 second point on first line
     * @param p3 first point on second line
     * @param p4 second point on second line
     * @param pInt1 point on first line that is closest to second line
     * @param s parameterized distance along first line from p1 to pInt1
     * @param pInt2 point on second line that is closest to first line
     * @param t parameterized distance along second line from p3 to pInt2
     * @return false if lines are parallel, true otherwise
     */
    static bool IntersectLines(
            const SimTK::Vec3& p1,
            const SimTK::Vec3& p2,
            const SimTK::Vec3& p3,
            const SimTK::Vec3& p4,
            SimTK::Vec3& pInt1,
            double& s,
            SimTK::Vec3& pInt2,
            double& t);

    /**
     * Compute the intersection of a line segment and a plane.
     *
     * @param pt1 first point on line
     * @param pt2 second point on line
     * @param plane normal vector of plane
     * @param d normal distance of plane to origin
     * @param inter intersection point of line and plane
     * @return true if line segment and plane intersect, false otherwise
     */
    static bool IntersectLineSegPlane(
            const SimTK::Vec3& pt1,
            const SimTK::Vec3& pt2,
            const SimTK::Vec3& plane,
            double d,
            SimTK::Vec3& inter);

    /**
     * Calculate the point (closestPt) on a line (linePt, line)
     * that is closest to a point (pt). 'line' does not need to
     * be normalized.
     *
     * @param pt the point
     * @param linePt a point on the line
     * @param line defines the line passing through linePt
     * @param closestPt the closest point
     * @param t parameterized distance from linePt along line to closestPt
     */
    static void GetClosestPointOnLineToPoint(
            const SimTK::Vec3& pt,
            const SimTK::Vec3& linePt,
            const SimTK::Vec3& line,
            SimTK::Vec3& closestPt,
            double& t);

    /**
     * Make a 3x3 direction cosine matrix for a rotation about the X axis.
     *
     * @param angle the rotation angle, in radians
     * @param m the 3x3 matrix
     */
    static void Make3x3DirCosMatrix(
            double angle,
            double mat[][3]);

    /**
     * Make a 4x4 direction cosine matrix from an axis/angle rotation.
     *
     * @param axis the axis of rotation
     * @param angle the angle, in radians
     * @param mat the matrix
     */
    static void ConvertAxisAngleTo4x4DirCosMatrix(
            const SimTK::Vec3& axis,
            double angle,
            double mat[][4]);

    /**
     * Compute the square of the distance between two points.
     *
     * @param point1 the first point
     * @param point2 the second point
     * @return the square of the distance
     */
    static double CalcDistanceSquaredBetweenPoints(
            const SimTK::Vec3& p1,
            const SimTK::Vec3& p2);

    /**
     * Compute the square of the distance between a point and a line.
     *
     * @param point the point
     * @param linePt a point on the line
     * @param line defines the line passing through linePt
     * @return the square of the distance
     */
    static double CalcDistanceSquaredPointToLine(
            const SimTK::Vec3& point,
            const SimTK::Vec3& linePt,
            const SimTK::Vec3& line);

    /**
     * Rotate a 4x4 transform matrix by 'angle' radians about axis 'axis'.
     *
     * @param matrix The 4x4 transform matrix
     * @param axis The axis about which to rotate
     * @param angle the amount to rotate, in radians
     */
    static void RotateMatrixAxisAngle(
            double matrix[][4],
            const SimTK::Vec3& axis,
            double angle);

    /**
     * Rotate a 4x4 transform matrix by 'angle' radians about axis 'axis'.
     *
     * TODO: legacy version: here for bootstrap testing
     *
     * @param matrix The 4x4 transform matrix
     * @param axis The axis about which to rotate
     * @param angle the amount to rotate, in radians
     */
    static void LegacyRotateMatrixAxisAngle(
            double matrix[][4],
            const SimTK::Vec3& axis,
            double angle);
};

/** @endcond **/

}; //namespace

#endif // __WrapMath_h__

