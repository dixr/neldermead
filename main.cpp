#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <limits>

constexpr const char* outfile_name = "neldermead.out";
constexpr const char* logfile_name = "neldermead.log";
constexpr const char* simplexfile_name = "neldermead.simplex";
constexpr const char* statefile_name = "neldermead.state"; // REFLECT, EXPAND, CONTRACT, SHRINK 1...N+1

constexpr double alpha = 1; // reflection (default are standard coefficients from wikipedia, also used in numerical recipes)
constexpr double gamma = 2; // expansion
constexpr double sigma = 0.5; // shrink
constexpr double rho = 0.5; // contraction
constexpr double eps = std::numeric_limits<double>::epsilon() * 1e6; // = TINY from numerical recipes

// not needed anymore:
// constexpr size_t max_iter = 1000; // max iterations for nelder mead
// constexpr bool accurate = true; // true: relative error criterion from numerical recipes, false: stopping criterion from brnt.eu (a very little tiny bit faster...)
//constexpr double epsilon = 1e-4; // meaning: if this is 1e-X and accurate=true, than all function values in the simplex will have the leading X significant digits in common! so 4 significant digits should be enough! (see optimize_gate/data/*)

using namespace std;

int main(int argc, char* argv[])
{
    // maintains neldermead.simplex, neldermead.out, neldermead.log
    // starts new optimization whenever there is no neldermead.simplex file (flushes neldermead.out)

    if(argc <= 2) {
        cerr << "usage: " << argv[0] << " parameters-file energy-file" << endl;
        return -1;
    }

    // read parameters (N)
    ifstream parametersfile(argv[1]);
    vector<double> parameters;
    if(!parametersfile) {
        cerr << "error: parameters-file '" << argv[1] << "' could not be opened" << endl;
        return -1;
    }
    double d;
    while(parametersfile >> p)
        parameters.push_back(p);
    const size_t N = parameters.size();
    if(N == 0){
       cerr << "error: parameters-file '" << argv[1] << "' has no data (only floating point numbers separated by whitespace are allowed)" << endl;
       return -1;
    }
    parametersfile.close();

    // read energy
    ifstream energyfile(argv[2]);
    double energy;
    if(!(energyfile >> energy)) {
        cerr << "error: energy-file '" << argv[1] << "' could not be read" << endl;
        return -1;
    }
    energyfile.close();

    // read simplex ? x (N+1) values (parameters, energy)
    ifstream simplexfilein(simplexfile_name);
    vector<double> simplex, fsimplex;
    for(int i = 0; simplexfilein >> d; ++i)
        if(i % (N+1) < N)
            simplex.push_back(p);
        else
            fsimplex.push_back(p);
    if(simplex.size() % N || simplex.size() > (N+1)*N || fsimplex.size() != N+1) {
        cerr << "error: simplex file '" << simplexfile_name << "' needs to contain ? x (N+1) values, i.e. ? times (parameters, energy), with ? <= N+1" << endl;
        return -1;
    }
    simplexfilein.close();

    // prepare output files (flushed if simplex was empty)
    ofstream outfile(outfile_name, simplex.empty() ? ios::out : ios::app);
    ofstream logfile(outfile_name, simplex.empty() ? ios::out : ios::app);
}

            size_t idx_min = 0, idx_max = 0, idx_second_max = 0;
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
                    output << "reached maximum number of iterations" << endl;
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
/**/
