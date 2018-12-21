#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <numeric>
#include <vector>
#include <limits>
#include <cmath>
#include <exception>

constexpr const char* logfile_name = "neldermead.log";
constexpr const char* statefile_name = "neldermead.state";

constexpr double alpha = 1;  // reflection (default coefficients from NR & wiki)
constexpr double eta = 2;  // expansion
constexpr double sigma = 0.5;  // shrink
constexpr double rho = 0.5;  // contraction
constexpr double eps = std::numeric_limits<double>::epsilon() * 1e6;  // = TINY
constexpr double inf = std::numeric_limits<double>::infinity();
constexpr double nanq = std::numeric_limits<double>::quiet_NaN();

/*  //program saves its current state in a state file of this format
    // (if not present, we start a new optimization, flushing the log file;
    // comments as shown here are not allowed in the file)
    EVALUATIONS 0
    state  // N+1 by N (x0; ...; xN)
    0 ... 0
    :     :
    0 ... 0
    FSIMPLEX  // N+1 (f(x0), ..., f(xN))
    0 ... 0
    X  // 3 by N (reflected, expanded, contracted)
    0 ... 0
    0 ... 0
    0 ... 0
    FX  // 3 (f(reflected), f(expanded), f(contracted))
    0 0 0
    MODE  // INITIALIZING idx, REFLECTING, EXPANDING, CONTRACTING, SHRINKING idx
*/

using namespace std;

void simplex_info(const vector<double>& simplex, const vector<double>& fsimplex,
                  const size_t N, size_t* idx_min, size_t* idx_max,
                  size_t* idx_second_max, double* delta,
                  vector<double>* centroid) {
    // partially order function values in fsimplex and compute delta
    vector<size_t> indices(N+1);
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(),
        [&fsimplex](size_t i, size_t j){return fsimplex[i] < fsimplex[j];});
    *idx_min = indices[0];
    *idx_max = indices[N];
    *idx_second_max = indices[N-1];
    *delta = 2*fabs(fsimplex[*idx_min] - fsimplex[*idx_max])
        / (fabs(fsimplex[*idx_min]) + fabs(fsimplex[*idx_max]) + eps);
        // relative error between simplex function values (from NR)

    // compute centroid of the simplex
    fill(centroid->begin(), centroid->end(), 0.);
    for (size_t i = 0; i < N+1; ++i)
        if (i != *idx_max) {
            for (size_t j = 0; j < N; ++j)
                (*centroid)[j] += simplex[i*N+j];
        }
    for (size_t j = 0; j < N; ++j)
        (*centroid)[j] /= N;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 3) {
            cerr << "usage: " << argv[0] << " parameters-file energy-file "
                 "[characteristic-scale-p0=1 ... characteristic-scale-pN-1=1]"
                 << endl;
            return -1;
        }

        // read parameters (also determines N)
        ifstream parametersfilein(argv[1]);
        vector<double> parameters;
        if (!parametersfilein) {
            cerr << "error: parameters-file '" << argv[1]
                 << "' could not be opened" << endl;
            return -1;
        }
        double d;
        while (parametersfilein >> d)
            parameters.push_back(d);
        parametersfilein.close();
        const size_t N = parameters.size();
        if (N < 2) {
           cerr << "error: parameters-file '" << argv[1]
                << "' has not enough data (only floating point numbers"
                " separated by whitespace are allowed; could only parse "
                << N << " parameters; need at least 2 parameters to optimize)"
                << endl;
           return -1;
        }

        // read energy
        ifstream energyfile(argv[2]);
        double energy;
        if (!(energyfile >> energy)) {
            cerr << "error: energy-file '" << argv[2]
                 << "' could not be read" << endl;
            return -1;
        }
        energyfile.close();

        // parse the characteristic scales for the parameters
        vector<double> parameter_scales(N, 1.);
        if (argc == 4) {
            parameter_scales.assign(N, stod(argv[3]));
        } else if (argc - 3U == N) {
            for (size_t i = 0; i < N; ++i)
                parameter_scales[i] = stod(argv[3+i]);
        } else if (argc != 3) {
            cerr << "error: wrong number of parameters given; either pass 1, N="
                 << N << ", or no characteristic scales" << endl;
            return -1;
        }

        // read state from statefile and prepare logfile
        ifstream statefilein(statefile_name);
        ofstream logfile;
        logfile.precision(10);
        vector<double> simplex((N+1)*N), fsimplex(N+1, energy);
        vector<double> x(3*N), fx(3, energy);  // reflected/expanded/contracted
        string mode = "INITIALIZING", part;
        size_t mode_idx = 0, evaluations = 1;
        if (statefilein >> part) {  // state file exists
            if (part != "EVALUATIONS" || !(statefilein >> evaluations))
                cerr << "error: state file '" << statefile_name
                     << "' is corrupt (part: EVALUATIONS)" << endl;
            ++evaluations;
            statefilein >> part;
            if (part != "SIMPLEX")
                cerr << "error: state file '" << statefile_name
                     << "' is corrupt (part: SIMPLEX)" << endl;
            for (double& d : simplex)
                statefilein >> d;
            statefilein >> part;
            if (part != "FSIMPLEX")
                cerr << "error: state file '" << statefile_name
                     << "' is corrupt (part: FSIMPLEX)" << endl;
            for (double& d : fsimplex)
                statefilein >> d;
            statefilein >> part;
            if (part != "X")
                cerr << "error: state file '" << statefile_name
                     << "' is corrupt (part: X)" << endl;
            for (double& d : x)
                statefilein >> d;
            statefilein >> part;
            if (part != "FX")
                cerr << "error: state file '" << statefile_name
                     << "' is corrupt (part: FX)" << endl;
            for (double& d : fx)
                statefilein >> d;
            statefilein >> mode;
            if (mode == "INITIALIZING" || mode == "SHRINKING") {
                if (!(statefilein >> mode_idx) || mode_idx > N)
                    cerr << "error: state file '" << statefile_name
                         << "' is corrupt: mode could not be parsed "
                         "or mode index " << mode_idx
                         << " is invalid (part: " << mode << ")" << endl;
            } else if (mode != "REFLECTING"
                    && mode != "EXPANDING"
                    && mode != "CONTRACTING") {
                cerr << "error: state file '" << statefile_name
                     << "' is corrupt: MODE '" << mode
                     << "' is unknown" << endl;
            }
            if (!statefilein) {
                cerr << "error: reached end of state file '" << statefile_name
                     << "'; it could not be parsed properly";
                return -1;
            }
            logfile.open(logfile_name, ios::app);
        } else {
            // state file did not exist, initialize simplex and x to the same
            // points so that each of them is valid; flush logfile too
            for (size_t i = 0; i < N+1; ++i)
                copy_n(parameters.begin(), N, simplex.begin() + i*N);
            for (size_t i = 0; i < 3; ++i)
                copy_n(parameters.begin(), N, x.begin() + i*N);
            logfile.open(logfile_name);
            logfile << "evaluation\tmode";
            for (size_t i = 0; i < N; ++i)
                logfile << "\tp" << i;
            logfile << "\tdelta\tenergy\n";
        }
        logfile << evaluations << '\t' << mode;
        for (double p : parameters)  // need to log parameters now;
            logfile << '\t' << p;    // they will be overwritten later
        statefilein.close();

        // main part of the algorithm
        size_t idx_min = 0, idx_max = 0, idx_second_max = 0;
        double delta = nanq;
        vector<double> centroid(N);
        if (mode == "INITIALIZING") {
            // building the simplex, currently at vector mode_idx
            copy_n(parameters.begin(), N, simplex.begin() + mode_idx*N);
            fsimplex[mode_idx] = energy;

            if (mode_idx < N) {
                copy_n(simplex.begin(), N, parameters.begin());
                parameters[mode_idx] += parameter_scales[mode_idx];
                ++mode_idx;
            } else {
                mode = "REFLECTING";
                simplex_info(simplex, fsimplex, N, &idx_min, &idx_max,
                             &idx_second_max, &delta, &centroid);
                for (size_t i = 0; i < N; ++i)
                    x[i] = centroid[i]
                         + alpha*(centroid[i] - simplex[idx_max*N+i]);
                copy_n(x.begin(), N, parameters.begin());
            }
        } else {
            simplex_info(simplex, fsimplex, N, &idx_min, &idx_max,
                         &idx_second_max, &delta, &centroid);
            // act on the other modes
            if (mode == "REFLECTING") {
                fx[0] = energy;
                if (fx[0] <= fsimplex[idx_min]) {
                    mode = "EXPANDING";
                    for (size_t i = 0; i < N; ++i)
                        x[N+i] = centroid[i] + eta*(x[i] - centroid[i]);
                    copy_n(x.begin() + N, N, parameters.begin());
                } else if (fx[0] >= fsimplex[idx_second_max]) {
                    mode = "CONTRACTING";
                    for (size_t i = 0; i < N; ++i)
                        x[2*N+i] = centroid[i]
                                 + rho*(simplex[idx_max*N+i] - centroid[i]);
                    copy_n(x.begin() + 2*N, N, parameters.begin());
                } else {
                    copy_n(parameters.begin(), N, simplex.begin() + idx_max*N);
                    fsimplex[idx_max] = energy;
                    mode = "REFLECTING";
                    simplex_info(simplex, fsimplex, N, &idx_min,
                                 &idx_max, &idx_second_max, &delta, &centroid);
                    for (size_t i = 0; i < N; ++i)
                        x[i] = centroid[i]
                             + alpha*(centroid[i] - simplex[idx_max*N+i]);
                    copy_n(x.begin(), N, parameters.begin());
                }
            } else if (mode == "EXPANDING") {
                fx[1] = energy;
                if (fx[0] <= fx[1]) {
                    copy_n(x.begin(), N, simplex.begin() + idx_max*N);
                    fsimplex[idx_max] = fx[0];
                } else {
                    copy_n(x.begin() + N, N, simplex.begin() + idx_max*N);
                    fsimplex[idx_max] = fx[1];
                }
                mode = "REFLECTING";
                simplex_info(simplex, fsimplex, N, &idx_min,
                             &idx_max, &idx_second_max, &delta, &centroid);
                for (size_t i = 0; i < N; ++i)
                    x[i] = centroid[i]
                         + alpha*(centroid[i] - simplex[idx_max*N+i]);
                copy_n(x.begin(), N, parameters.begin());
            } else if (mode == "CONTRACTING") {
                fx[2] = energy;
                if (fx[2] <= fsimplex[idx_max]) {
                    copy_n(x.begin() + 2*N, N, simplex.begin() + idx_max*N);
                    fsimplex[idx_max] = fx[2];
                    mode = "REFLECTING";
                    simplex_info(simplex, fsimplex, N, &idx_min, &idx_max,
                                 &idx_second_max, &delta, &centroid);
                    for (size_t i = 0; i < N; ++i)
                        x[i] = centroid[i]
                             + alpha*(centroid[i] - simplex[idx_max*N+i]);
                    copy_n(x.begin(), N, parameters.begin());
                } else {
                    mode = "SHRINKING";
                    mode_idx = (idx_min > 0 ? 0 : 1);
                    for (size_t i = 0; i < N+1; ++i)
                        if (i != idx_min) {
                            for (size_t j = 0; j < N; ++j)
                                simplex[i*N+j] = simplex[idx_min*N+j]
                                               + sigma*(simplex[i*N+j]
                                                      - simplex[idx_min*N+j]);
                        }
                    copy_n(simplex.begin(), N, parameters.begin());
                }
            } else if (mode == "SHRINKING") {
                fsimplex[mode_idx] = energy;
                if (++mode_idx == idx_min)
                    ++mode_idx;
                if (mode_idx < N) {
                    copy_n(simplex.begin() + mode_idx*N, N, parameters.begin());
                } else {
                    mode = "REFLECTING";
                    simplex_info(simplex, fsimplex, N, &idx_min, &idx_max,
                                 &idx_second_max, &delta, &centroid);
                    for (size_t i = 0; i < N; ++i)
                        x[i] = centroid[i]
                             + alpha*(centroid[i] - simplex[idx_max*N+i]);
                    copy_n(x.begin(), N, parameters.begin());
                }
            }
        }

        // update and close logfile
        logfile << '\t' << delta << '\t' << energy << '\n';
        logfile.close();

        // save new parameters
        ofstream parametersfileout(argv[1]);
        parametersfileout.precision(10);
        for (size_t i = 0; i < N; ++i)
            parametersfileout << parameters[i] << (i % 2 == 1 ? '\n' : ' ');
        parametersfileout.close();

        // save state file
        ofstream statefileout(statefile_name);
        statefileout.precision(10);
        statefileout << "EVALUATIONS " << evaluations << '\n';
        statefileout << "SIMPLEX\n";
        for (size_t i = 0; i < simplex.size(); ++i)
            statefileout << simplex[i] << (i % N == N-1 ? '\n' : '\t');
        statefileout << "FSIMPLEX\n";
        for (size_t i = 0; i < fsimplex.size(); ++i)
            statefileout << fsimplex[i] << (i == N ? '\n' : '\t');
        statefileout << "X\n";
        for (size_t i = 0; i < x.size(); ++i)
            statefileout << x[i] << (i % N == N-1 ? '\n' : '\t');
        statefileout << "FX\n";
        for (size_t i = 0; i < fx.size(); ++i)
            statefileout << fx[i] << (i == 2 ? '\n' : '\t');
        statefileout << mode;
        if (mode == "INITIALIZING" || mode == "SHRINKING")
            statefileout << ' ' << mode_idx;
        statefileout << '\n';
        statefileout.close();

        // save the best point so far as minimum
        // (start with the simplex as it always contains valid points)
        simplex.insert(simplex.end(), x.begin(), x.end());
        fsimplex.insert(fsimplex.end(), fx.begin(), fx.end());
        size_t argmin = min_element(fsimplex.begin(), fsimplex.end())
                      - fsimplex.begin();
        ofstream(string("min-")+argv[2]) << setprecision(10)
                                         << fsimplex[argmin] << endl;
        ofstream minparametersfileout(string("min-")+argv[1]);
        minparametersfileout.precision(10);
        for (size_t i = 0; i < N; ++i)
            minparametersfileout << simplex[argmin*N+i]
                                 << (i % 2 == 1 ? '\n' : ' ');
        minparametersfileout.close();
    }
    catch(const std::invalid_argument& e) {
        cerr << "error: invalid argument (" << e.what()
             << "): some of the characteristic "
             "scales are not convertible to double" << endl;
        return -1;
    }
    catch(const std::exception& e) {
        cerr << "error (exception): " << e.what() << endl;
        return -1;
    }
}
