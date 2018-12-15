#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <numeric>
#include <vector>
#include <limits>
#include <cstdlib>
#include <exception>

constexpr const char* logfile_name = "neldermead.log";
constexpr const char* statefile_name = "neldermead.state";

constexpr double alpha = 1; // reflection (default are standard coefficients from Numerical Recipes, also on wikipedia)
constexpr double gamma = 2; // expansion
constexpr double sigma = 0.5; // shrink
constexpr double rho = 0.5; // contraction
constexpr double eps = std::numeric_limits<double>::epsilon() * 1e6; // = TINY from numerical recipes
constexpr double inf = std::numeric_limits<double>::infinity();
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
constexpr double pi = 3.14159'26535'89793'23846'26433'83279'502884;

/* // program saves its current state in a state file of this format
   // (if not present, we start a new optimization, flushing the log file; comments as shown here are not allowed in the file)
    EVALUATION 0
    state // N+1 by N (x0; ...; xN)
    0 ... 0
    :     :
    0 ... 0
    FSIMPLEX // N+1 (f(x0), ..., f(xN))
    0 ... 0
    X // 3 by N (reflected, expanded, contracted)
    0 ... 0
    0 ... 0
    0 ... 0
    FX // 3 (f(reflected), f(expanded), f(contracted))
    0 0 0
    MODE // can be INITIALIZING idx, REFLECTING, EXPANDING, CONTRACTING, SHRINKING idx
*/

using namespace std;

void simplex_info(const vector<double>& simplex, const vector<double> fsimplex, const size_t N,
                  size_t& idx_min, size_t idx_max, size_t& idx_second_max, double& delta, vector<double>& centroid)
{
    // partially order function in fsimplex values and compute delta
    vector<size_t> indices(N+1);
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&fsimplex](size_t i, size_t j){return fsimplex[i] < fsimplex[j];});
    idx_min = indices[0];
    idx_max = indices[N];
    idx_second_max = indices[N-1];
    delta = 2*abs(fsimplex[idx_min] - fsimplex[idx_max])/(abs(fsimplex[idx_min]) + abs(fsimplex[idx_max]) + eps);

    // compute centroid of the simplex
    fill(centroid.begin(), centroid.end(), 0.);
    for(size_t i = 0; i < N+1; ++i)
        if(i != idx_max)
            for(size_t j = 0; j < N; ++j)
                centroid[j] += simplex[i*N+j];
    for(size_t j = 0; j < N; ++j)
        centroid[j] /= N;
}

int main(int argc, char* argv[])
{
    try
    {
        if(argc <= 2) {
            cerr << "usage: " << argv[0] << " parameters-file energy-file [characteristic-parameter-scale=PI]" << endl;
            return -1;
        }

        // if no scale is given, use pi to initialize the simplex; change that for other parameters!
        const double parameter_scale = (argc == 3 ? stod(argv[2]) : pi);

        // read parameters (N)
        ifstream parametersfilein(argv[1]);
        vector<double> parameters;
        if(!parametersfilein) {
            cerr << "error: parameters-file '" << argv[1] << "' could not be opened" << endl;
            return -1;
        }
        double d;
        while(parametersfilein >> d)
            parameters.push_back(d);
        const size_t N = parameters.size();
        if(N < 2){
           cerr << "error: parameters-file '" << argv[1] << "' has not enough data (only floating point numbers separated by whitespace are allowed; could only parse " << N << " parameters; need at least 2 parameters to optimize)" << endl;
           return -1;
        }
        parametersfilein.close();

        // read energy
        ifstream energyfile(argv[2]);
        double energy;
        if(!(energyfile >> energy)) {
            cerr << "error: energy-file '" << argv[1] << "' could not be read" << endl;
            return -1;
        }
        energyfile.close();

        // read state from statefile and prepare logfile
        ifstream statefilein(statefile_name);
        ofstream logfile;
        vector<double> simplex((N+1)*N), fsimplex(N+1);
        vector<double> x(3*N), fx(3); // reflected, expanded, contracted
        string mode = "INITIALIZING", part;
        size_t mode_idx = 0, evaluation = 1;
        if(statefilein >> part) // state file exists
        {
            if(part != "EVALUATION" || !(statefilein >> evaluation)) cerr << "error: state file '" << statefile_name << "' is corrupt (part: EVALUATION)" << endl;
            statefilein >> part;
            if(part != "SIMPLEX") cerr << "error: state file '" << statefile_name << "' is corrupt (part: SIMPLEX)" << endl;
            for(double& d : simplex)
                statefilein >> d;
            statefilein >> part;
            if(part != "FSIMPLEX") cerr << "error: state file '" << statefile_name << "' is corrupt (part: FSIMPLEX)" << endl;
            for(double& d : fsimplex)
                statefilein >> d;
            statefilein >> part;
            if(part != "X") cerr << "error: state file '" << statefile_name << "' is corrupt (part: X)" << endl;
            for(double& d : x)
                statefilein >> d;
            statefilein >> part;
            if(part != "FX") cerr << "error: state file '" << statefile_name << "' is corrupt (part: FX)" << endl;
            for(double& d : fx)
                statefilein >> d;
            statefilein >> mode;
            if(mode == "INITIALIZING" || mode == "SHRINKING") {
                if(!(statefilein >> mode_idx) || mode_idx > N)
                    cerr << "error: state file '" << statefile_name << "' is corrupt: mode could not be parsed or mode index " << mode_idx << " is invalid (part: " << mode << ")" << endl;
            }
            else if(mode != "REFLECTING" && mode != "EXPANDING" && mode != "CONTRACTING")
                cerr << "error: state file '" << statefile_name << "' is corrupt: MODE '" << mode << "' is unknown" << endl;
            if(!statefilein) {
                cerr << "error: reached end of state file '" << statefile_name << "'; it could not be parsed properly";
                return -1;
            }
            logfile.open(logfile_name, ios::app);
        }
        else // state file did not exist; flush logfile too
        {
            logfile.open(logfile_name);
            logfile << "evaluation\tmode";
            for(size_t i = 0; i < N; ++i)
                logfile << "\tp" << i;
            logfile << "\tdelta\tenergy" << endl;
        }
        logfile << evaluation << '\t' << mode;
        for(double p : parameters) // need to log parameters now; they will be overwritten later
            logfile << '\t' << p;
        statefilein.close();

        // main part of the algorithm
        size_t idx_min = 0, idx_max = 0, idx_second_max = 0;
        double delta = nan;
        vector<double> centroid(N);
        if(mode == "INITIALIZING") // building the simplex, currently at vector mode_idx
        {
            copy_n(parameters.begin(), N, simplex.begin() + mode_idx*(N+1));
            fsimplex[mode_idx] = energy;

            if(mode_idx < N)
            {
                copy_n(simplex.begin(), N, parameters.begin());
                parameters[mode_idx] += parameter_scale; // NOTE: maybe should do fmod?
                ++mode_idx;
            }
            else
            {
                mode = "REFLECTING";
                simplex_info(simplex, fsimplex, N, idx_min, idx_max, idx_second_max, delta, centroid);
                for(size_t i = 0; i < N; ++i)
                    x[i] = centroid[i] + alpha*(centroid[i] - simplex[idx_max*N+i]);
                copy_n(x.begin(), N, parameters.begin());
            }
        }
        else
        {
            simplex_info(simplex, fsimplex, N, idx_min, idx_max, idx_second_max, delta, centroid);

            // act on mode
            if(mode == "REFLECTING")
            {
                fx[0] = energy;
                if(fx[0] <= fsimplex[idx_min])
                {
                    mode = "EXPANDING";
                    for(size_t i = 0; i < N; ++i)
                        x[N+i] = centroid[i] + gamma*(x[i] - centroid[i]);
                    copy_n(x.begin() + N, N, parameters.begin());
                }
                else if(fx[0] >= fsimplex[idx_second_max])
                {
                    mode = "CONTRACTING";
                    for(size_t i = 0; i < N; ++i)
                        x[2*N+i] = centroid[i] + rho*(simplex[idx_max*N+i] - centroid[i]);
                    copy_n(x.begin() + 2*N, N, parameters.begin());
                }
                else
                {
                    copy_n(parameters.begin(), N, simplex.begin() + idx_max*(N+1));
                    fsimplex[idx_max] = energy;

                    mode = "REFLECTING";
                    simplex_info(simplex, fsimplex, N, idx_min, idx_max, idx_second_max, delta, centroid);
                    for(size_t i = 0; i < N; ++i)
                        x[i] = centroid[i] + alpha*(centroid[i] - simplex[idx_max*N+i]);
                    copy_n(x.begin(), N, parameters.begin());
                }
            }
            else if(mode == "EXPANDING")
            {
                fx[1] = energy;
                if(fx[0] <= fx[1])
                {
                    copy_n(x.begin(), N, simplex.begin() + idx_max*(N+1));
                    fsimplex[idx_max] = fx[0];
                }
                else
                {
                    copy_n(x.begin() + N, N, simplex.begin() + idx_max*(N+1));
                    fsimplex[idx_max] = fx[1];
                }

                mode = "REFLECTING";
                simplex_info(simplex, fsimplex, N, idx_min, idx_max, idx_second_max, delta, centroid);
                for(size_t i = 0; i < N; ++i)
                    x[i] = centroid[i] + alpha*(centroid[i] - simplex[idx_max*N+i]);
                copy_n(x.begin(), N, parameters.begin());
            }
            else if(mode == "CONTRACTING")
            {
                fx[2] = energy;
                if(fx[2] <= fsimplex[idx_max])
                {
                    copy_n(x.begin() + 2*N, N, simplex.begin() + idx_max*(N+1));
                    fsimplex[idx_max] = fx[2];

                    mode = "REFLECTING";
                    simplex_info(simplex, fsimplex, N, idx_min, idx_max, idx_second_max, delta, centroid);
                    for(size_t i = 0; i < N; ++i)
                        x[i] = centroid[i] + alpha*(centroid[i] - simplex[idx_max*N+i]);
                    copy_n(x.begin(), N, parameters.begin());
                }
                else
                {
                    mode = "SHRINKING";
                    mode_idx = (idx_min > 0 ? 0 : 1);

                    for(size_t i = 0; i < N+1; ++i)
                        if(i != idx_min)
                            for(size_t j = 0; j < N; ++j)
                                simplex[i*N+j] = simplex[idx_min*N+j] + sigma*(simplex[i*N+j] - simplex[idx_min*N+j]);
                    copy_n(simplex.begin(), N, parameters.begin());
                }
            }
            else if(mode == "SHRINKING")
            {
                fsimplex[mode_idx] = energy;

                if(++mode_idx == idx_min)
                    ++mode_idx;

                if(mode_idx < N)
                    copy_n(simplex.begin() + mode_idx*N, N, parameters.begin());
                else
                {
                    mode = "REFLECTING";
                    simplex_info(simplex, fsimplex, N, idx_min, idx_max, idx_second_max, delta, centroid);
                    for(size_t i = 0; i < N; ++i)
                        x[i] = centroid[i] + alpha*(centroid[i] - simplex[idx_max*N+i]);
                    copy_n(x.begin(), N, parameters.begin());
                }
            }
        }

        // update and close logfile
        logfile << '\t' << delta << '\t' << energy << endl;
        logfile.close();

        // save new parameters
        ofstream parametersfileout(argv[1]);
        for(size_t i = 0; i < N; ++i)
            parametersfileout << parameters[i] << (i % 2 == 1 ? '\n' : ' ');
        parametersfileout.close();

        // save state file
        ofstream statefileout(statefile_name);
        statefileout << "EVALUATION " << ++evaluation << endl;
        statefileout << "SIMPLEX" << endl;
        for(size_t i = 0; i < simplex.size(); ++i)
            statefileout << simplex[i] << (i % N == N-1 ? '\n' : '\t');
        statefileout << "FSIMPLEX" << endl;
        for(size_t i = 0; i < fsimplex.size(); ++i)
            statefileout << fsimplex[i] << (i == N ? '\n' : '\t');
        statefileout << "X" << endl;
        for(size_t i = 0; i < x.size(); ++i)
            statefileout << x[i] << (i % N == N-1 ? '\n' : '\t');
        statefileout << "FX" << endl;
        for(size_t i = 0; i < fx.size(); ++i)
            statefileout << fx[i] << (i == 2 ? '\n' : '\t');
        statefileout << mode;
        if(mode == "INITIALIZING" || mode == "SHRINKING")
            statefileout << ' ' << mode_idx;
        statefileout << endl;
        statefileout.close();
    }
    catch(const std::exception& e)
    {
        cerr << "error (exception): " << e.what() << std::endl;
        return -1;
    }
}

/*
// constexpr size_t max_iter = 1000; // max iterations for nelder mead
// constexpr bool accurate = true; // true: relative error criterion from numerical recipes, false: stopping criterion from brnt.eu (a very little tiny bit faster...)
//constexpr double epsilon = 1e-4; // meaning: if this is 1e-X and accurate=true, than all function values in the simplex will have the leading X significant digits in common! so 4 significant digits should be enough! (see optimize_gate/data/)

            size_t num_eval = 0;
            size_t num_iter = 0;
            double delta = numeric_limits<double>::infinity();

            ofstream output(outfilename, ios::app);
            output.precision(10);

            Eigen::Array<double,N+1,1> fsimplex; // contains f(x0); ...; f(xN)
            Eigen::Array<double,3,N> x; // contains x_reflection; x_expansion; x_contraction
            Eigen::Array<double,3,1> fx; // contains f(x_reflection); f(x_expansion); f(x_contraction)
            enum : size_t { refl = 0, expd, cont }; // indices for reflected, expanded and contracted point

            for(size_t i = 0; i < N+1; ++i) // removed parallel region here and in shrink since objective() should decide on its own which of its steps need to be parallelized (note: omp_get_thread_num() return 0 in nested regions was no error since the innermost parallel team had indeed only one thread; results were also fine)
                fsimplex[i] = objective(simplex.row(i));
            num_eval += N+1;

            while(true)
            {
                if(!simplexfilename.empty())
                {
                    ofstream simplexfile(simplexfilename); simplexfile.precision(17);
                    for(size_t i = 0; i < N+1; ++i)
                        for(size_t j = 0; j < N; ++j)
                            simplexfile << simplex(i,j) << (j == N-1 ? '\n' : '\t');
                }

                output << "[iter=" << num_iter << ":eval=" << num_eval << "]" << endl;
                for(size_t i = 0; i < N+1; ++i)
                    output << "x" << i << " = (" << simplex.row(i) << "),\tf(x" << i << ") = " << fsimplex[i] << endl;

                // partially order function num_eval
                array<size_t,N+1> indices;
                iota(indices.begin(), indices.end(), 0);
                sort(indices.begin(), indices.end(), [&fsimplex](size_t i, size_t j){return fsimplex[i] < fsimplex[j];});
                idx_min = indices[0];
                idx_max = indices[N];
                idx_second_max = indices[N-1];
                output << "found idx_min=" << idx_min << ", idx_max=" << idx_max << ", idx_second_max=" << idx_second_max << endl;

                // check for convergence
                if(accurate)
                    delta = 2*fabs(fsimplex[idx_min] - fsimplex[idx_max])/(fabs(fsimplex[idx_min]) + fabs(fsimplex[idx_max]) + eps); // stopping criterion from numerical recipes
                else
                    delta = sqrt((fsimplex - fsimplex.sum()/(N+1)).square().sum()/(N+1)); // stopping criterion from brnt.eu
                output << "delta = " << delta << endl;
                if(delta < epsilon)
                {
                    output << "converged!" << endl;
                    break;
                }
                else if(++num_iter >= max_iter)
                {
                    output << "reached maximum number of evaluations" << endl;
                    break;
                }

                Eigen::Array<double,1,N> centroid = Eigen::Array<double,1,N>::Zero();
                for(size_t i = 0; i < N+1; ++i)
                    if(i != idx_max)
                        centroid += simplex.row(i);
                centroid /= N;
                output << "centroid = (" << centroid << ")" << endl;

                // try reflection
                x.row(refl) = centroid + alpha*(centroid - simplex.row(idx_max));
                fx[refl] = objective(x.row(refl));
                ++num_eval;
                output << "reflected = (" << x.row(refl) << "),\tf(reflected) = " << fx[refl] << endl;

                if(fx[refl] <= fsimplex[idx_min])
                {
                    // if reflected is so good, try expanding a bit more
                    x.row(expd) = centroid + gamma*(x.row(refl) - centroid);
                    fx[expd] = objective(x.row(expd));
                    ++num_eval;
                    output << "expanded = (" << x.row(expd) << "),\tf(expanded) = " << fx[expd] << endl;

                    if(fx[expd] < fx[refl])
                    {
                        output << "==> expanded is best, so EXPAND\n" << endl;
                        simplex.row(idx_max) = x.row(expd);
                        fsimplex[idx_max] = fx[expd];
                    }
                    else
                    {
                        output << "==> reflected point is best, so REFLECT\n" << endl;
                        simplex.row(idx_max) = x.row(refl);
                        fsimplex[idx_max] = fx[refl];
                    }
                }
                else if(fx[refl] >= fsimplex[idx_second_max])
                {
                    // reflection brings no improvement, try contracting
                    x.row(cont) = centroid + rho*(simplex.row(idx_max) - centroid);
                    fx[cont] = objective(x.row(cont));
                    ++num_eval;
                    output << "contracted = (" << x.row(cont) << "),\tf(contracted) = " << fx[cont] << endl;

                    if(fx[cont] < fsimplex[idx_max])
                    {
                        output << "==> contraction is better than reflection, so CONTRACT\n" << endl;
                        simplex.row(idx_max) = x.row(cont);
                        fsimplex[idx_max] = fx[cont];
                    }
                    else
                    {
                        output << "==> all new points are bad, so SHRINK\n" << endl;
                        for(size_t i = 0; i < N+1; ++i)
                            if(i != idx_min)
                            {
                                simplex.row(i) = simplex.row(idx_min) + sigma*(simplex.row(i) - simplex.row(idx_min));
                                fsimplex[i] = objective(simplex.row(i));
                            }
                        num_eval += N;
                    }
                }
                else
                {
                    output << "==> reflected point is ok, so REFLECT\n" << endl;
                    simplex.row(idx_max) = x.row(refl);
                    fsimplex[idx_max] = fx[refl];
                }
            }

            output << "==> final solution is x = (" << simplex.row(idx_min) << ") with f(x) = " << fsimplex[idx_min] << endl;
            if(pfresult) *pfresult = fsimplex[idx_min];
            if(pdelta) *pdelta = delta;
            if(pnum_iter) *pnum_iter = num_iter;
            if(pnum_eval) *pnum_eval = num_eval;
            return simplex.row(idx_min);
        }

        template<int N, class... U>
        Eigen::Array<double,1,N> run_initdel(
            const Eigen::Array<double,1,N>& init, // contains initial guess to build simplex
            const Eigen::Array<double,1,N>& del, // characteristic length scales in each dimension (used to construct initial simplex as recommended by numerical recipes)
            U&&... u)
        {
            Eigen::Array<double,N+1,N> simplex = init.matrix().replicate(N+1,1); // containts x0; ...; xN
            for(size_t i = 1; i < N+1; ++i)
                simplex.row(i)[i-1] += del[i-1];
            return run_simplex<N>(simplex, forward<U>(u)...);
        }

        template<int N, typename objective_type, class... U>
        Eigen::Array<double,1,N> resume_run(
            const objective_type& objective,
            const string& simplexfilename = "",
            U&&... u)
        {
            Eigen::Array<double,N+1,N> simplex;
            ifstream simplexfile(simplexfilename);
            for(size_t i = 0; i < N+1; ++i)
                for(size_t j = 0; j < N; ++j)
                    if(!(simplexfile >> simplex(i,j)))
                        throw runtime_error("[optimization::nelder_mead::resume_run] Could not load simplex from file '"+simplexfilename+"'; error in row "+to_string(i)+" column "+to_string(j));
            return run_simplex<N>(simplex, objective, simplexfilename, forward<U>(u)...);
        }
    }
}
*/
