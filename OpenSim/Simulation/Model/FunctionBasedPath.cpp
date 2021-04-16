#include "FunctionBasedPath.h"

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PointBasedPath.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>

#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>

// Returns `true` if changing the supplied `Coordinate` changes the moment arm
// of the supplied `PointBasedPath` (PBP)
static bool coordinateAffectsPBP(
    OpenSim::PointBasedPath const& pbp,
    OpenSim::Coordinate const& c,
    SimTK::State& state) {

    static constexpr int numCoordProbingSteps = 3;
    static constexpr double minMomentArmChangeRequired = 0.001;

    bool initialLockedState = c.getLocked(state);
    double initialValue = c.getValue(state);

    c.setLocked(state, false);

    double start = c.getRangeMin();
    double end = c.getRangeMax();
    double step = (end - start) / numCoordProbingSteps;

    bool affectsCoord = false;
    for (double v = start; v <= end; v += step) {
        c.setValue(state, v);
        double ma = pbp.computeMomentArm(state, c);

        if (std::abs(ma) >= minMomentArmChangeRequired) {
            affectsCoord = true;
            break;
        }
    }

    c.setLocked(state, initialLockedState);
    c.setValue(state, initialValue);

    return affectsCoord;
}

// Fills the output vector with n linearly-spaced values in the inclusive range
// [begin, end]
static void linspace(double begin,
                     double end,
                     size_t n,
                     std::vector<double>& out) {
    if (n == 0) {
        return;
    }

    if (n == 1) {
        out.push_back(begin);
        return;
    }

    double step = (end - begin) / (n-1);
    for (size_t i = 0; i < n-1; ++i) {
        out.push_back(begin + i*step);
    }
    out.push_back(end);
}

// A discretization of n points along the (incluside) interval [begin, end]
struct Discretization final {
    double begin;
    double end;
    int nPoints;
    double gridsize;
};

// overload for linspace that is specialized for `Discretization`s
static void linspace(const Discretization& d, std::vector<double>& out) {
    assert(d.nPoints >= 0);

    out.clear();
    linspace(d.begin, d.end, static_cast<size_t>(d.nPoints), out);
}

class Interpolate final {
private:
    // the dimensionality of the interpolation (i.e. a 2D fit == 2)
    int dimensions;

    // vector containing the orthogonal discretization vectors
    // e.g. disc[0] = xrange, disc[1] = yrange etc.
    std::vector<std::vector<double>> discretizations;

    std::vector<double> discSizes;

    // array containing the evaluations over the discretizations
    std::vector<double> evals;

    // vector containing index closest to discretization point
    std::vector<int> n;

    // vector containing fraction of closest index to discretization
    // point to next
    std::vector<double> u;

    // array of polynomial evaluation
    std::vector<std::vector<double>> beta;

    // vector of an index in the evalsPair
    std::vector<int> loc;

    // OPENSIM INTEGRATION
    std::vector<Discretization> dS;
    std::vector<const OpenSim::Coordinate *> coords;

public:
    Interpolate() = default;

    // Precomputed data constructor
    Interpolate(
            std::vector<OpenSim::Coordinate const*> coordsIn,
            std::vector<Discretization> dSIn,
            std::vector<double> evalsIn) :

        dimensions(static_cast<int>(dSIn.size())),
        evals(evalsIn),
        n(dimensions, 0),
        u(dimensions, 0),
        loc(dimensions, 0),
        dS(dSIn),
        coords(coordsIn)
    {
        std::vector<double> dc_;
        for (size_t i = 0; i < static_cast<size_t>(dimensions); i++) {
            beta.push_back({0, 0, 0, 0});

            linspace(dS[i], dc_);
            discretizations.push_back(dc_);
            discSizes.push_back(dS[i].gridsize);
        }
    }

    // Constructor for non FBP/PBP related interpolation
    Interpolate(
            std::vector<std::vector<double>> discretizationIn,
            std::vector<std::pair<std::vector<int>,double>> evalsPair) :

        dimensions(static_cast<int>(discretizationIn.size())),
        discretizations(discretizationIn),
        n(dimensions,0),
        u(dimensions,0),
        loc(dimensions,0)
    {
        assert(discretizations.size() == evalsPair[0].first.size());

        // allow it to work with the new struct method
        for (int i = 0; i < dimensions; i++) {
            Discretization dc;
            dc.begin = discretizations[i][0];
            dc.end = discretizations[i][discretizations[i].size()-1];
            dc.nPoints = static_cast<int>(discretizations[i].size());
            dc.gridsize = (dc.end - dc.begin)/(dc.nPoints-1);

            dS.push_back(dc);
        }

        for (int i = 0; i < dimensions; i++) {
            beta.push_back({0,0,0,0});
            discSizes.push_back(discretizations[i].size());
        }

        // I'm aware this is still stupid but it is what it is for now
        int numOfLoops = 1;
        for (int i = 0; i < dimensions; i++) {
            numOfLoops *= discretizations[i].size();
        }

        std::vector<int> cnt(dimensions, 0);
        for (int i = 0; i < numOfLoops; i++) {
            for (unsigned k = 0; k < evalsPair.size(); k++) {
                if (evalsPair[k].first == cnt) {
                    evals.push_back(evalsPair[k].second);
                }
            }

            for (int x = dimensions - 1; x >= 0; x--) {
                if (cnt[x] != dS[x].nPoints-1) {
                    cnt[x] += 1;
                    break;
                }

                if (cnt[x] == dS[x].nPoints-1) {
                    for (int y = x; y < dimensions; y++) {
                        cnt[y] = 0;
                    }
                }
            }
        }
    }

    // General interface constructor
    // allows one to create an interface constructor as shown below
    // with a vector of coordinates
    Interpolate(OpenSim::PointBasedPath const& pbp,
                OpenSim::Coordinate const** cBegin,
                OpenSim::Coordinate const** cEnd,
                SimTK::State& st,
                std::vector<int>& nPoints) :

        dimensions(static_cast<int>(nPoints.size())),
        n(dimensions,0),
        u(dimensions,0),
        loc(dimensions,0)
    {
        // get the size of the 'vector'. Inclusive bottom, exclusive top
        std::ptrdiff_t n = (cEnd - cBegin) + 1;
        assert(n > 0);
        assert(dimensions == (int)n);

        // put all coordinate pointers in a vector to later unpack an incoming
        // state to a vector of coordinate values
        for (int i = 0; i < dimensions; i++) {
            coords.push_back(cBegin[i]);
        }

        // unlock coordinates
        for (int i = 0; i < dimensions; i++) {
            const OpenSim::Coordinate& c = *cBegin[i];
            c.setLocked(st, false);
        }

        // make discretization objects for interpolation class instance
        Discretization dc;
        for (int i = 0; i < dimensions; i++) {
            const OpenSim::Coordinate& c = *cBegin[i];
            dc.begin = std::max(c.getRangeMin(), -static_cast<double>(SimTK_PI));
            dc.end = std::min(c.getRangeMax(), static_cast<double>(SimTK_PI));
            dc.nPoints = nPoints[i];
            dc.gridsize = (dc.end - dc.begin) / (dc.nPoints - 1);
            dS.push_back(dc);
        }

        // slightly extend the bound for accurate interpolation on the edges
        for (int dim = 0; dim < dimensions; dim++) {
            dS[dim].begin   -= 1*dS[dim].gridsize;
            dS[dim].end     += 2*dS[dim].gridsize;
            dS[dim].nPoints += 3;
        }

        // compute the total number of evaluations we need to do
        int numOfLoops = 1;
        for (int i = 0; i < dimensions; i++) {
            numOfLoops *= dS[i].nPoints;
        }

        std::vector<int> cnt(dimensions);
        std::vector<double> coordValues(dimensions);
        for (int i = 0; i < numOfLoops; i++) {
            for (int ii = 0; ii < dimensions; ii++) {
                coordValues[ii] = dS[ii].begin + (cnt[ii]*dS[ii].gridsize);
                const OpenSim::Coordinate& c = *cBegin[ii];
                c.setValue(st,coordValues[ii]);
            }
            evals.push_back((pbp.getLength(st)));

            // update cnt values
            for (int x = dimensions-1; x >= 0; x--) {
                if (cnt[x] != dS[x].nPoints-1){
                    cnt[x] += 1;
                    break;
                } else {
                    for (int y = x; y < dimensions; y++) {
                        cnt[y] = 0;
                    }
                }
            }
        }

        // just make it for using the old getInterp method
        std::vector<double> dc_;
        for (int i = 0; i < dimensions; i++) {
            beta.push_back({0,0,0,0});

            dc_.clear();
            linspace(dS[i].begin, dS[i].end, dS[i].nPoints, dc_);
            discretizations.push_back(dc_);
            discSizes.push_back(dS[i].nPoints);
        }
    }

    // Basic constructor having a vector of coordinate pointers as input
    Interpolate(
            OpenSim::PointBasedPath const& pbp,
            std::vector<OpenSim::Coordinate const*> coords,
            SimTK::State& st,
            std::vector<int>& nPoints)  : Interpolate(pbp,
                                                      &coords[0],
                                                      &coords[coords.size()-1],
                                                      st,
                                                      nPoints)
    {
    }

    std::vector<double> getRange(int i) const {
        return discretizations[i];
    }

    int getDimension() const {
        return dimensions;
    }

    std::vector<Discretization> getdS() const {
        return dS;
    }

    std::vector<double> getEvals() const {
        return evals;
    }

    double getEval() {
        // get the wrapping length evaluation given a vector 'loc' which contains,
        // in ascending dimension, the index in each dimension
        int factor = 1;
        int idx = 0;

        for (int i = 0; i < dimensions-1; i++) {
            factor = 1;
            for (int ii = i+1; ii <= dimensions-1; ii++) {
                factor *= discSizes[ii];
            }
            idx += loc[i]*factor;
        }
        idx += loc[loc.size()-1];

        assert(idx <= evals.size());

        return evals[idx];
    }

    // 0th derivative
    double getInterp(const std::vector<double>& x) {
        // This is the main interpolation function
        // IN:  x, a vector of points within the considered interpolation range
        // OUT: eval, the interpolated value
        assert(x.size() == dimensions);

        // get the index of the closest range value to the discretization point
        for (int i = 0; i < dimensions; i++) {
            auto it = std::find_if(std::begin(discretizations[i]),
                                   std::end(discretizations[i]),
                                   [&](double j){return j >= x[i];});

            n[i] = static_cast<int>(std::distance(discretizations[i].begin(), it) - 1);
        }

        // compute remaining fraction
        for (int i = 0; i < dimensions; i++) {
            u[i] = (x[i]-discretizations[i][n[i]])/
                    (discretizations[i][2]-discretizations[i][1]);
        }

        // compute the polynomials (already evaluated)
        for (int i = 0; i < dimensions; i++) {
            // compute binomial coefficient
            beta[i][0] = 0.5*(u[i]-1)*(u[i]-1)*(u[i]-1)*u[i]*(2*u[i]+1);
            beta[i][1] = -0.5*(u[i] - 1)*(6*u[i]*u[i]*u[i]*u[i] -
                                          9*u[i]*u[i]*u[i] + 2*u[i] + 2);
            beta[i][2] = 0.5*u[i]*(6*u[i]*u[i]*u[i]*u[i] - 15*u[i]*u[i]*u[i] +
                                   9*u[i]*u[i] + u[i] + 1);
            beta[i][3] = -0.5*(u[i] - 1)*u[i]*u[i]*u[i]*(2*u[i] - 3);
        }

        // loop over all the considered points (n-dimensional) and multiply the
        // evaluation with the weight

        // create an array containing the  lengths of each discretization direction
        int discrLoopCnt[dimensions];
        for (int i = 0; i < dimensions; i++) {
            discrLoopCnt[i] = -1;
        }

        double z = 0;
        bool breakWhile = false;
        bool allTrue;
        double Beta = 1;
        while (discrLoopCnt[0] < 3) {
            Beta = 1;
            for (int i = 0; i < dimensions; i++) {
                Beta = Beta * beta[i][discrLoopCnt[i]+1];
            }

            for (int i = 0; i < dimensions; i++) {
                loc[i] = discrLoopCnt[i] + n[i];
            }

            z += getEval() * Beta;

            // from the back to the front, check if we're already at the maximum
            // iteration on that 'nested' for loop or else increment with 1.
            // In short, everything starts with [-1,-1,-1,...] and we keep adding
            // ones until the array of the loops becomes [ 2, 2, 2, ...]
            for (int x = dimensions-1; x >= 0; x--) {
                if (discrLoopCnt[x] != 2) {
                    discrLoopCnt[x] += 1;
                    break;
                }
                if (discrLoopCnt[x] == 2) {
                    for (int y = x; y < dimensions; y++) {
                        discrLoopCnt[y] = -1;
                    }
                }
            }

            // checking exit conditions
            if (breakWhile) {
                break;
            }

            // loop through to check whether all are at max cnt or not
            allTrue = true;
            for (int i = 0; i < dimensions; i++) {
                if (discrLoopCnt[i] != 2){
                    allTrue = false;
                    break;
                }
            }
            // if all true (all are 2) set breakWhile to break on the next iteration
            if (allTrue) {
                breakWhile = true;
            }
        }
        return z;
    }

    // getLength with mapping from State to vector x
    double getInterp(const SimTK::State& s) {
        assert(coords.size() != 0);

        std::vector<double> x;
        for (const OpenSim::Coordinate* c : coords) {
            x.push_back(c->getValue(s));
        }

        return getInterp(x);
    }

    // getLength based on the state
    double getLength(const SimTK::State& s) {
        return getInterp(s);
    }

    // 1st derivative
    double getInterpDer(const std::vector<double>& x,
                        int coordinate,
                        double h = 0.0001) {

        assert(x.size() == dimensions);
        assert(coordinate <= dimensions-1);
        assert(h>0 && h < (dS[coordinate].end - x[coordinate]));

        // This is the main interpolation function for derivatives
        // IN:  x, a vector of points within the considered interpolation range
        //      coordinate, the generalized coordinate of which we take derivative
        // OUT: eval, the interpolated value
        double f1 = getInterp(x);

        std::vector<double> x2 = x;
        x2[coordinate] += h;
        double f2 = getInterp(x2);

        return (f2 - f1)/h;
    }

    double getInterpDer(const SimTK::State& s, int coordinate) {
        assert(coords.size() != 0);

        std::vector<double> x;
        for (size_t i = 0; i < coords.size(); i++) {
            x.push_back(coords[i]->getValue(s));
        }
        return getInterpDer(x,coordinate);
    }

    double getInterpDer(const SimTK::State& s,
                        const OpenSim::Coordinate& aCoord) {

        for (size_t i = 0; i < coords.size(); i++) {
            if (coords[i]->getName() == aCoord.getName()) {
                return getInterpDer(s, i);
            }
        }

        return 0.0;
    }

    double getLengtheningSpeed(const SimTK::State& s) {
        double lengtheningSpeed = 0.0;

        for (int i = 0; i < dimensions; i++) {
            double firstDeriv = getInterpDer(s, i);
            double coordinateSpeedVal = coords[i]->getSpeedValue(s);

            lengtheningSpeed += firstDeriv * coordinateSpeedVal;
        }

        return lengtheningSpeed;
    }
};

// ----- FunctionBasedPath implementation -----

struct OpenSim::FunctionBasedPath::Impl final {
    mutable Interpolate interp;

    Impl* clone() const {
        auto p = std::unique_ptr<Impl>{new Impl{}};
        p->interp = interp;
        return p.release();
    }
};

OpenSim::FunctionBasedPath OpenSim::FunctionBasedPath::fromPointBasedPath(
        const Model& model,
        const PointBasedPath& pbp)
{
    FunctionBasedPath fbp;

    // copy relevant data from source PBP
    fbp.upd_Appearance() = pbp.get_Appearance();
    fbp.setPathPointSet(pbp.getPathPointSet());
    fbp.setPathWrapSet(pbp.getWrapSet());

    std::unique_ptr<OpenSim::Model> modelClone{model.clone()};
    SimTK::State& initialState = modelClone->initSystem();
    modelClone->equilibrateMuscles(initialState);
    modelClone->realizeVelocity(initialState);

    std::vector<const Coordinate *> coordsThatAffectPBP;
    for (Coordinate const& c : modelClone->getComponentList<Coordinate>()){
        if (c.getMotionType() == Coordinate::MotionType::Coupled) {
            continue;
        }

        if (!coordinateAffectsPBP(pbp, c, initialState)) {
            continue;
        }

        coordsThatAffectPBP.push_back(&c);
        std::cout << "affecting coord: " << c.getName() << std::endl;
    }

    // create vector for number of interpolation points
    std::vector<int> nPoints(coordsThatAffectPBP.size(), 5);

    // reinitialize the default-initialized interp with the points
    fbp._impl->interp = Interpolate{
            pbp,
            std::move(coordsThatAffectPBP),
            initialState,
            nPoints
    };

    return fbp;
}

static Interpolate readInterp(const std::string &path) {
    // data file for FunctionBasedPath top-level structure:
    //
    //     - plaintext
    //     - line-oriented
    //     - 3 sections (HEADERS, DISCRETIZATIONS, and EVALUATIONS) separated
    //       by blank lines
    //
    // details:
    //
    // line                                  content
    // -------------------------------------------------------------------
    // 0                                     ID (string)
    // 1                                     <empty>
    // 2                                     DIMENSIONS (int)
    // 3                                     <empty>
    // --- end HEADERS
    // 4                                     DISCRETIZATION_1
    // 5                                     DISCRETIZATION_2
    // 6-                                    DISCRETIZATION_n
    // 3+DIMENSIONS                          DISCRETIZATION_${DIMENSIONS}
    // 3+DIMENSIONS+1                        <empty>
    // --- end DISCRETIZATIONS
    // 3+DIMENSIONS+2                        EVAL_1_1
    // 3+DIMENSIONS+3                        EVAL_1_2
    // 3+DIMENSIONS+4-                       EVAL_1_m
    // 3+DIMENSIONS+${points_in_range}       EVAL_1_${points_in_range}
    // 3+DIMENSIONS+${points_in_range}+1     EVAL_2_1
    // ...
    // ? (depends on points_in_range for each discretization)
    // ?                                     EVAL_${DIMENSIONS}_${points_in_range}
    // EOF
    //
    // where DISCRETIZATION_$n is tab-delimited (\t) line that describes
    // evenly-spaced points along the `n`th DIMENSION:
    //
    // column           content
    // ------------------------------------------------
    // 0                range_start (float)
    // 1                range_end (float)
    // 2                points_in_range (int)
    // 3                grid_size (float)
    //
    // and EVAL_$n_$m is a line that contains a single real-valued evaluation
    // of the curve at:
    //
    //     range_start + (m * ((range_end - range_start)/points_in_range))
    //
    // along DIMENSION `n`

    std::ifstream file{path};

    if (!file) {
        std::stringstream msg;
        msg << path << ": error opening FunctionBasedPath data file";
        OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
    }

    int lineno = 0;
    std::string line;

    // ID (ignored)
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()) {
            break;
        }
    }

    // DIMENSIONS
    int dimensions = -1;
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()) {
            if (dimensions <= -1) {
                std::stringstream msg;
                msg << path << ": L" << lineno << ": unexpected blank line (expected dimensions)";
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }
            break;
        }

        dimensions = std::stoi(line.c_str());
    }

    // [DISCRETIZATION_1..DISCRETIZATION_n]
    std::vector<Discretization> discretizations;
    discretizations.reserve(static_cast<size_t>(dimensions));
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()) {
            if (discretizations.size() != static_cast<size_t>(dimensions)) {
                std::stringstream msg;
                msg << path << ": L" << lineno << ": invalid number of discretizations in this file: expected: " << dimensions << " got " << discretizations.size();
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }

            break;  // end of DISCRETIZATIONs
        }

        // DISCRETIZATION_i
        Discretization d;
        std::sscanf(line.c_str(), "%lf\t%lf\t%i\t%lf", &d.begin, &d.end, &d.nPoints, &d.gridsize);
        discretizations.push_back(d);
    }

    // [EVAL_1_1..EVAL_n_m]
    std::vector<double> evals;
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()){
            break;  // end of EVALs
        }

        evals.push_back(atof(line.c_str()));
    }

    // sanity check
    {
        size_t expectedEvals = 0;
        for (const Discretization& d : discretizations) {
            expectedEvals += static_cast<size_t>(d.nPoints);
        }

        if (expectedEvals != evals.size()) {
            std::stringstream msg;
            msg << path << ": L" << lineno << ": invalid number of function evaluations in this file: expected " << expectedEvals << " got " << evals.size();
            OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
        }
    }

    // TODO: figure out which coordinates affect this path - must be serialized
    //       in the file

    return Interpolate{
        std::vector<const OpenSim::Coordinate*>(static_cast<size_t>(dimensions)),
        std::move(discretizations),
        std::move(evals)
    };
}

OpenSim::FunctionBasedPath OpenSim::FunctionBasedPath::fromDataFile(const std::string &path)
{
    FunctionBasedPath fbp;
    fbp.setDataPath(path);

    // the data file effectively only contains data for Interpolate, data file
    // spec is in there
    fbp._impl->interp = readInterp(path);

    return fbp;
}

OpenSim::FunctionBasedPath::FunctionBasedPath() : GeometryPath{}, _impl{new Impl{}}
{
    constructProperty_data_path("");
}
OpenSim::FunctionBasedPath::FunctionBasedPath(const FunctionBasedPath&) = default;
OpenSim::FunctionBasedPath::FunctionBasedPath(FunctionBasedPath&&) = default;
OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(const FunctionBasedPath&) = default;
OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(FunctionBasedPath&&) = default;
OpenSim::FunctionBasedPath::~FunctionBasedPath() noexcept = default;

void OpenSim::FunctionBasedPath::extendFinalizeFromProperties() {
    bool hasBackingFile = !get_data_path().empty();
    bool hasInMemoryFittingData =  !_impl->interp.getdS().empty();

    if (hasBackingFile) {
        // load the file, always
        //
        // TODO: this should ideally be cached
        _impl->interp = readInterp(get_data_path());
    } else if (hasInMemoryFittingData) {
        // do nothing: just use the already-loaded in-memory fitting data
        return;
    } else {
        // error: has no backing file and has no in-memory fitting data
        std::stringstream msg;
        msg << getAbsolutePathString() << ": cannot call `.finalizeFromProperties()` on this `" << getConcreteClassName() << "`: the path has no fitting associated with it. A `FunctionBasedPath` must either have a valid `data_path` property that points to its underlying fitting data data, or be initialized using `FunctionBasedPath::fromPointBasedPath` (i.e. generate new fitting data)";
        OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
    }
}

double OpenSim::FunctionBasedPath::getLength(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _lengthCV)) {
        return getCacheVariableValue(s, _lengthCV);
    }

    // else: compute it
    double rv = _impl->interp.getLength(s);
    setCacheVariableValue(s, _lengthCV, rv);
    return rv;
}

void OpenSim::FunctionBasedPath::setLength(const SimTK::State& s, double length) const
{
    setCacheVariableValue(s, _lengthCV, length);
}

double OpenSim::FunctionBasedPath::getLengtheningSpeed(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _speedCV)) {
        return getCacheVariableValue(s, _speedCV);
    }

    // else: compute it
    double rv = _impl->interp.getLengtheningSpeed(s);
    setCacheVariableValue(s, _speedCV, rv);
    return rv;
}

void OpenSim::FunctionBasedPath::setLengtheningSpeed(const SimTK::State& s, double speed) const
{
    setCacheVariableValue(s, _speedCV, speed);
}

double OpenSim::FunctionBasedPath::computeMomentArm(const SimTK::State& s,
                                                    const Coordinate& aCoord) const
{
    return _impl->interp.getInterpDer(s,aCoord);
}

void OpenSim::FunctionBasedPath::printContent(std::ostream& out) const
{
    // see FunctionBasedPath::fromDataFile(const std::string &path) for format
    // spec

    // HEADERS
    out << get_data_path() << '\n'
        << '\n'
        << _impl->interp.getDimension() << '\n'
        << '\n';

    // DISCRETIZATIONS
    for (const Discretization& d : _impl->interp.getdS()) {
        out << d.begin << '\t'
            << d.end << '\t'
            << d.nPoints << '\t'
            << d.gridsize << '\n';
    }
    out << '\n';

    // EVALS
    for (double d : _impl->interp.getEvals()) {
        out << d << '\n';
    }

    out.flush();
}