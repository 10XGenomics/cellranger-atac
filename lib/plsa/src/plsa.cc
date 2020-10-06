// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
#include "docopt.h"
#include "cmdargs.h"

int OMP_THREADS = std::min(4, omp_get_num_procs());

// sparse matrix class
template <class T>
class sparse_matrix{
    public:
        std::vector<int64_t> rowidx, colidx;
        std::vector<T> vals_csr, vals_csc;
        std::vector<int> cols, rows;
        int nr;
        int nc;
        std::vector<int64_t> perm;

    // constructor using precomputed rows
    sparse_matrix(){}

    // read csv file
    sparse_matrix(const std::string filename){
        std::string extension = get_extension(filename);
        if (extension == "csv")
            read_dense_csv_into_csr(filename, vals_csr, rowidx, cols, nc, nr);
        else if (extension == "mtx")
            read_mtx_into_csr(filename, vals_csr, rowidx, cols, nc, nr);
        else
            ExitWithMsg("Error opening file, unidentified file extension", extension);

        // convert to csc format, and store the permutation from csr to csc
        // eats up some memory
        std::vector<int> rrows(vals_csr.size(),0);
        for(int r = 0; r < nr; r++)
            for(int64_t x = rowidx[r]; x < rowidx[r+1]; x++)
                rrows[x] = r;
        std::vector<std::vector<int64_t> > colsvals(nc);
        for(int64_t x = 0; x < vals_csr.size(); x++)
            colsvals[cols[x]].push_back(x);

        perm.reserve(vals_csr.size());
        rows.reserve(vals_csr.size());
        vals_csc.reserve(vals_csr.size());
        colidx.resize(nc+1,0);
        for(int c = 0; c < nc; c++){
            perm.insert( perm.end(), colsvals[c].begin(), colsvals[c].end() );
            for (auto x : colsvals[c]){
                rows.push_back(rrows[x]);
                vals_csc.push_back(vals_csr[x]);
            }
            colidx[c+1] = colsvals[c].size() + colidx[c];
        }
        std::cout << "Matrix dimensions: " << nr << " x " << nc << std::endl;
        std::vector<std::vector<int64_t> >().swap(colsvals);
    }

    // destructor
    ~sparse_matrix(){
        std::vector<int64_t>().swap(rowidx);
        std::vector<int64_t>().swap(colidx);
        std::vector<int>().swap(rows);
        std::vector<int>().swap(cols);
        std::vector<T>().swap(vals_csr);
        std::vector<T>().swap(vals_csc);
        std::vector<int64_t>().swap(perm);
        nr = 0, nc = 0;
    }

    int ncols(){ return nc;}
    int nrows(){ return nr;}
    int64_t nnz() {return vals_csr.size();}

    private:

    std::string get_extension(const std::string s){
       size_t i = s.rfind('.', s.length());
       if (i != std::string::npos) {
          return(s.substr(i+1, s.length() - i));
       }
       return("");
    }

    void read_dense_csv_into_csr(const std::string filename, std::vector<T>& valscsr,
                                 std::vector<int64_t>& rowid, std::vector<int>& colscsr, 
                                 int& NC, int& NR){
        // assumes rows are barcodes, columns are peaks
        std::cout << "Reading csv assuming one row per sample...";
        std::ifstream file(filename);
        std::string line, val;

        // read in csr format
        rowid.push_back(0);
        while (std::getline(file, line))
        {
            std::stringstream iss(line);
            int col = 0;
            while(std::getline(iss, val, ',')){
                double v = std::stof(val);
                if (v > 1e-28){
                    valscsr.push_back(v);
                    colscsr.push_back(col);
                }
                col++;
            }
            NC = col;
            rowid.push_back(valscsr.size());
        }
        NR = rowid.size()-1;
        std::cout << "Done!" << std::endl;
    }

    void read_mtx_into_csr(const std::string filename, std::vector<T>& valscsr,
                                 std::vector<int64_t>& rowid, std::vector<int>& colscsr, 
                                 int& NC, int& NR){
        // first column is 1-index into peaks
        // second column is 1-index into barcodes
        // third column is the value
        std::cout << "Reading mtx assuming {column 1: col index, column 2: row index, column 3: val}...";
        std::ifstream file(filename);
        std::string line, val;

        // read in coo format, get it zero indexed
        struct entry{
            int r; int c; T v;
        };
        std::vector<entry> entries;
        int line_num = 0;
        while (std::getline(file, line))
        {
            line_num++;
            if(line_num < 3) continue;
            std::stringstream iss(line);
            int colm = 0;
            while(std::getline(iss, val, ' ')){
                double v = std::stof(val);
                if(line_num == 3){
                    if (colm == 0)
                        NC = v;
                    if (colm == 1)
                        NR = v;
                    if (colm == 2)
                        entries.resize(v);
                }else{
                    if (colm == 0)
                        entries[line_num-4].c = (v-1);
                    if (colm == 1)
                        entries[line_num-4].r = (v-1);
                    if (colm == 2)
                        entries[line_num-4].v = v;
                }
                colm++;
            }
        }
        // sort entries by 'r'
        std::sort(entries.begin(), entries.end(), [](const entry & e1, const entry & e2){ if(e1.r == e2.r) return (e1.c < e2.c); return (e1.r < e2.r);});
        // convert coo to csr format
        valscsr.reserve(entries.size());
        colscsr.reserve(entries.size());
        for(const auto& e: entries){
            valscsr.push_back(e.v);
            colscsr.push_back(e.c);
        }
        rowid.reserve(NR+1);
        int curr_r = 0;
        for(int64_t i = 0; i < entries.size(); i++){
            if (entries[i].r == curr_r){
                rowid.push_back(i);
                curr_r++;
            }
        }
        rowid.push_back(entries.size());

        std::cout << "Done!" << std::endl;
    }
};

// dense matrix class in row-major
#define ROWMAJOR true
#define COLMAJOR false
typedef bool format;
template <class T>
class dense_matrix{
    private:
        std::vector<T> vals;
        int nr;
        int nc;
        format major;

    public:
    // initializers
    dense_matrix(){}

    dense_matrix(int r, int c, T def, format M){
        nr = r;
        nc = c;
        vals.resize(r * c, def);
        major = M;
    }

    dense_matrix(int r, int c, format M){
        nr = r;
        nc = c;
        vals.resize(r * c, 0);
        major = M;
    }

    // accessors
    T get(int r, int c) const {
        if (major == ROWMAJOR)
            return vals[r * nc + c];
        else
            return vals[c * nr + r];
    }

    int nrows() const{
        return nr;
    }

    int ncols() const{
        return nc;
    }

    // modifiers
    void set(int r, int c, T v){
        if (major == ROWMAJOR)
            vals[r * nc + c] = v;
        else
            vals[c * nr + r] = v;
    }

    // copy
    dense_matrix<T> (const dense_matrix<T>& source){
        if (this == &source)
            return;
        nr = source.nr;
        nc = source.nc;
        major = source.major;
        vals = source.vals;
    }

    // assignment
    dense_matrix<T>& operator= (const dense_matrix<T>& source){
        if(this == &source)
            return *this;
        nr = source.nr;
        nc = source.nc;
        major = source.major;
        vals = source.vals;
        return *this;
    } 

    // destuctor
    ~dense_matrix(){
        nr = 0;
        nc = 0;
        std::vector<T>().swap(vals);
    }
};

// pLSA class
template <class T>
class pLSA{

     private:
         int D, W, K;
         sparse_matrix<T> N_dw;
         std::vector<T> N_d;
         dense_matrix<double> P_ZD, P_WZ;
         std::vector<double> P_Z;

     public:

         pLSA(const std::string matrix_file, int ntopics=30){
             // get sparse out
             N_dw = sparse_matrix<T>(matrix_file);
             K = ntopics;
             D = N_dw.nrows();
             W = N_dw.ncols();
             // populate N_d from N_dw
             N_d.resize(D,0);
             for(int d = 0; d < D; d++)
                 for(int64_t ix = N_dw.rowidx[d]; ix < N_dw.rowidx[d+1]; ix++)
                     N_d[d] += N_dw.vals_csr[ix];
         }

         // perform iterations
         void fit_transform( int maxiter=100, double tol=0.001){

             int niter = (maxiter > 1) ? maxiter : 1;
             double tolerance = tol;

             // initialize P_wz, P_zd
             dense_matrix<double> P_wz(K, W, 0, ROWMAJOR);
             omp_set_dynamic(0);
             omp_set_num_threads(OMP_THREADS);
             #pragma omp parallel for
             for(int z = 0; z < K; z++){
                 double sum = 0;
                 for(int w = 0; w < W; w++){
                     P_wz.set(z,w, z + w + 1); //randint(100) + 1;
                     sum += P_wz.get(z,w);
                 }
                 for(int w = 0; w < W; w++){
                     P_wz.set(z,w, P_wz.get(z,w) / sum);
                 }
             }
             dense_matrix<double> P_zd(D, K, 0, ROWMAJOR);
             omp_set_dynamic(0);
             omp_set_num_threads(OMP_THREADS);
             #pragma omp parallel for
             for(int d = 0; d < D; d++){
                 double sum = 0;
                 for(int z = 0; z < K; z++){
                     P_zd.set(d,z, d + z + 1); //randint(100) + 1;
                     sum += P_zd.get(d,z);
                 }
                 for(int z = 0; z < K; z++){
                     P_zd.set(d,z, P_zd.get(d,z) / sum);
                 }
             }

             // iterate
             IterateEM(tolerance, niter, P_wz, P_zd);

             // save data
             double N = 0.0;
             for(auto nd : N_d)
                 N += nd;

             P_Z.resize(K, 0);
             P_ZD = P_zd;
             omp_set_dynamic(0);
             omp_set_num_threads(OMP_THREADS);
             #pragma omp parallel for
             for(int z = 0; z < K; z++){
                 for(int d = 0; d < D; d++){
                    P_ZD.set(d,z, P_ZD.get(d,z) * N_d[d] / N);
                    P_Z[z] += P_ZD.get(d,z);
                 }
             }
             P_WZ = P_wz;
         }

         void Write(std::string& path){
             // write out P(d,z)
             std::ofstream file(path + "/transformed_matrix.csv");
             for(int r = 0; r < D; r++){
                 for(int c = 0; c < K; c++){
                     file << std::setprecision(12) << P_ZD.get(r,c);
                     if(c<K-1) file << ",";
                 }
                 file << "\n";
             }
             file.close();

             // write out P(w|z)
             std::ofstream file2(path + "/components.csv");
             for(int r = 0; r < K; r++){
                 for(int c = 0; c < W; c++){
                     file2 << std::setprecision(12) << P_WZ.get(r,c);
                     if(c<W-1) file2 << ",";
                 }
                 file2 << "\n";
             }
             file2.close();

             // write out P(z)
             std::ofstream file3(path + "/topic_relevance.csv");
             for(int r = 0; r < K; r++){
                 file3 << P_Z[r];
                 if (r < K-1)
                     file3 << std::setprecision(12) << std::endl;
             }
             file3.close();
         }


    private:
         void IterateEM(double tol, int niter,
                 dense_matrix<double>& P_wz,
                 dense_matrix<double>& P_zd){

             // track tolerance and reduce Expectation step
             dense_matrix<double> P_zd_old = P_zd;
             dense_matrix<double> P_wz_old = P_wz;
             int64_t nnz = N_dw.nnz();

             // allocate P_zdw_sum
             std::vector<double> P_zdw_sum_csr(N_dw.nnz(), 0);
             std::vector<double> P_zdw_sum_csc(N_dw.nnz(), 0);

             std::cout << "================================" << std::endl;
             std::cout << "Number of iterations: " << niter << std::endl;
             std::cout << "Number of topics: " << K << std::endl;
             std::cout << "Sparsity factor: " << nnz * 1.0 / (D*W) << std::endl;
             std::cout << "Number of cells: " << D << std::endl;
             std::cout << "Number of sites: " << W << std::endl;
             std::cout << "================================" << std::endl;

             double res = 10000;
             for(int iter = 0; iter < niter; iter++){
                 std::cout << "Running iteration: " << iter << std::endl;

                 // Expectation step
                 omp_set_dynamic(0);
                 omp_set_num_threads(OMP_THREADS);
                 #pragma omp parallel for
                 for(int w = 0; w < W; w++){
                     for(int64_t iy = N_dw.colidx[w]; iy < N_dw.colidx[w+1]; iy++){
                         P_zdw_sum_csc[iy] = 0;
                         for(int z = 0; z < K; z++)
                             P_zdw_sum_csc[iy] += P_zd_old.get(N_dw.rows[iy],z) * P_wz_old.get(z,w);
                         P_zdw_sum_csr[N_dw.perm[iy]] = P_zdw_sum_csc[iy];
                     }
                 }

                 // Maximization step
                 omp_set_dynamic(0);
                 omp_set_num_threads(OMP_THREADS);
                 #pragma omp parallel for
                 for(int w = 0; w < W; w++){
                     for(int z = 0; z < K; z++){
                         P_wz.set(z,w, 0);
                         for(int64_t iy = N_dw.colidx[w]; iy < N_dw.colidx[w+1]; iy++)
                             P_wz.set(z,w, P_wz.get(z,w) + 
                                             N_dw.vals_csc[iy] * P_zd_old.get(N_dw.rows[iy],z) * P_wz_old.get(z,w) / P_zdw_sum_csc[iy]);
                     }
                 }
                 omp_set_dynamic(0);
                 omp_set_num_threads(OMP_THREADS);
                 #pragma omp parallel for
                 for(int z = 0; z < K; z++){
                     double sum = 0;
                     for(int w = 0; w < W; w++)
                         sum += P_wz.get(z,w);
                     for(int w = 0; w < W; w++)
                         P_wz.set(z,w, P_wz.get(z,w) / sum);
                 }

                 omp_set_dynamic(0);
                 omp_set_num_threads(OMP_THREADS);
                 #pragma omp parallel for
                 for(int d = 0; d < D; d++){
                     for(int z = 0; z < K; z++){
                         P_zd.set(d,z, 0);
                         for(int64_t ix = N_dw.rowidx[d]; ix < N_dw.rowidx[d+1]; ix++)
                             P_zd.set(d,z, P_zd.get(d,z) + N_dw.vals_csr[ix] * P_zd_old.get(d,z) * P_wz_old.get(z,N_dw.cols[ix]) / P_zdw_sum_csr[ix]);
                         P_zd.set(d,z, P_zd.get(d,z) / N_d[d]);
                     }
                 }

                 res = matrix_norm(P_zd, P_zd_old);
                 std::cout << "Residue: " << res << std::endl;
                 if (res < tol && iter > 10 ){
                     std::cout << "Converged!" << std::endl;
                     break;
                 }
                 // copy can be avoided by referencing alternately
                 P_zd_old = P_zd;
                 P_wz_old = P_wz;
             }
         }

         // the 1-norm
         double matrix_norm(dense_matrix<double>& A, dense_matrix<double>& B){
             std::vector<double> rowsums(A.nrows(),0);
             omp_set_dynamic(0);
             omp_set_num_threads(OMP_THREADS);
             #pragma omp parallel for
             for(int r = 0; r < A.nrows(); r++)
                 for(int c = 0; c < A.ncols(); c++)
                     rowsums[r] += fabs(A.get(r,c) - B.get(r,c));
             double max = 0;
             for(auto c: rowsums)
                 max = std::max(max,c);
             return max;
         }
};

void process_args(argmap& args,  command_args& cmdline){
     auto& MATRIX_FILE = cmdline.MATRIX_FILE;
     auto& OUT = cmdline.OUT;
     auto& TOPICS = cmdline.TOPICS;
     auto& ITER = cmdline.ITER;
     auto& TOL = cmdline.TOL;
     auto& THREADS = cmdline.THREADS;

     // process inputs
     typedef typename dc_tag<int>::type itag;
     typedef typename dc_tag<double>::type dtag;
     typedef typename dc_tag<std::string>::type stag;
     process_arg(MATRIX_FILE, args, "<infile>", stag());
     process_arg(OUT, args, "<outdir>", stag());
     process_arg(TOPICS, args, "--topics", itag());
     process_arg(ITER, args, "--iter", itag());
     process_arg(TOL, args, "--tol", dtag());
     process_arg(THREADS, args, "--nt", itag());

     // sanity check
     if (!file_exists(MATRIX_FILE)) ExitWithMsg("matrix file not found:", MATRIX_FILE);
     if (!directory_exists(OUT)) ExitWithMsg("output directory does not exist:", OUT);
     if (TOPICS < 10){
         std::cout << "Defaulting to 10 hidden topics" << std::endl;
         TOPICS = 10;
     }
     if (ITER < 10){
         std::cout << "Minimum iterations = 10" << std::endl;
         ITER = 10;
     }
     if(TOL < 0){
         std::cout << "tolerance cannot be negative. Defaulting to 3E-3" << std::endl;
         TOL = 3E-3;
     }
     if (THREADS > omp_get_num_procs() || THREADS < 1){
         std::cout << "thread request for " << omp_get_num_procs() 
                   << " threads is invalid. Defaulting to min(4, available procs)" << std::endl;
         THREADS = std::min(4, omp_get_num_procs());
     }
     cmdline.print();
}

static const char USAGE[] =
R"(plsa

    Usage:
      plsa <infile> <outdir> [--topics=<t>] [--iter=<i>] [--tol=<o>] [--nt=<n>]
      plsa (-h | --help)
      plsa --version

    Options:
      -h --help     Show this screen.
      --version     Show version.
      --topics=<t>  Number of hidden topics [default: 10]
      --iter=<i>    Max iters [default: 2500]
      --tol=<o>     Tolerance for convergence [default: 0.003]
      --nt=<n>      Number of threads [default: 4]
)";

int main( int argc, char *argv[] )
{
     // process args
     std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
                         { argv + 1, argv + argc },
                         true,               // show help if requested
                         "Probabilistic Latent Semantic Analysis 2.0");  // version string
     command_args cmdline;
     process_args(args, cmdline);
     OMP_THREADS = cmdline.THREADS;

     double start = omp_get_wtime();

     // initialize plsa
     pLSA<double> plsa(cmdline.MATRIX_FILE, cmdline.TOPICS);

     // Do iterations
     plsa.fit_transform(cmdline.ITER, cmdline.TOL);

     // Write out
     plsa.Write(cmdline.OUT);

     std::cout << "Time elapsed: " << omp_get_wtime() - start << " seconds" << std::endl;
}
