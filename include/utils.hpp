#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <filesystem>
#include <cfloat>

inline std::string construct_metadata_path( std::string primary_identifier, std::string secondary_identifier, const int m, const int n, const int k, const int s, const int seed ){
    std::string metadata_path = primary_identifier;
    metadata_path += "/" + secondary_identifier;
    metadata_path += "/m" + std::to_string(m) + "_n" + std::to_string(n) + "_k" + std::to_string(k);
    metadata_path += "/scl" + std::to_string(s);
    metadata_path += "/seed" + std::to_string(seed);
    return metadata_path;
}

template <typename T>
void calculate_relative_errors( const int m, const int n, const double *C, const int ldc1, const T *C_ref, const int ldc2, double *rel_errors, const int ldc3 ){
    const long double flush_threshold = static_cast<long double>(DBL_TRUE_MIN);
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            const long double ci  = static_cast<long double>(C[i + j*ldc1]);
            const long double cri = static_cast<long double>(C_ref[i + j*ldc2]);
            const long double abs_err = std::abs(ci - cri);
            const long double denom   = std::abs(ci) + std::abs(cri);
            const long double rel_err_ld = abs_err == 0.0L ? 0.0L : abs_err / denom;
            const double rel_err = static_cast<double>(
                rel_err_ld != 0.0L && rel_err_ld < flush_threshold
                ? flush_threshold
                : rel_err_ld
            );
            if ( rel_err > rel_errors[i + j*ldc3] ){
                rel_errors[i + j*ldc3] = rel_err;
            }
        }
    }
}

template <typename T>
void calculate_f_of_n( const int m, const int n, const double *C, const int ldc1, const T *C_ref, const int ldc2, const T *C_ref_abs, const int ldc3, double *f_of_n, const int ldc4 ){
    
    const long double unit_roundoff = std::ldexp(1.0L, -53);
    
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            const long double ci    = static_cast<long double>(C[i + j*ldc1]);
            const long double cri   = static_cast<long double>(C_ref[i + j*ldc2]);
            const long double cri_a = static_cast<long double>(C_ref_abs[i + j*ldc3]);
            const long double abs_err = std::abs(ci - cri);
            f_of_n[i + j*ldc4] = static_cast<double>(abs_err == 0.0L ? 0.0L : abs_err / (cri_a * unit_roundoff));
        }
    }
}

template <typename T>
void dump_matrix(const int m, const int n, const T *A, const int lda, std::string outfile_path ){

    std::filesystem::create_directories(std::filesystem::path(outfile_path).parent_path());

    std::ofstream binfile(outfile_path, std::ios::binary);
    for ( int j = 0; j < n; ++j ){
        binfile.write(reinterpret_cast<const char*>(A + j * lda), m * sizeof(T));
    }
    binfile.close();
}

template <typename T>
int load_matrix( std::string infile_path, const int m, const int n, T *A, const int lda ){

    std::ifstream binfile(infile_path, std::ios::binary);
    if ( !binfile.good() ){
        return 0;
    }

    for ( int j = 0; j < n; ++j ){
        binfile.read(reinterpret_cast<char*>(A + j * lda), m * sizeof(T));
        if ( binfile.gcount() != static_cast<std::streamsize>(m * sizeof(T)) ){
            fprintf(stderr,
                "ERROR: short read from '%s' at column %d (expected %zu bytes, got %zu)\n",
                infile_path.c_str(), j, m * sizeof(T), (size_t)binfile.gcount());
            exit(1);
        }
    }
    binfile.close();

    return 1;
}

template <typename T_to, typename T_from>
void copy_matrix( const int m, const int n, T_to *A, const int lda, const T_from *B, const int ldb ){
    
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[i + j*lda] = B[i + j*ldb];
        }
    }

    return;
}

template <typename T>
void abs_matrix( const int m, const int n, T *A, const int lda ){
    
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[i + j*lda] = std::abs(A[i + j*lda]);
        }
    }

    return;
}

template <typename T>
void permute_rows( const int seed, const int m, const int n, T *A, const int lda ){
    std::vector<int> p(m);
    for ( int i = 0; i < m; ++i ){
        p[i] = i;
    }
    std::mt19937 gen(seed);
    std::shuffle(p.begin(), p.end(), gen);

    std::vector<T> A_tmp(m * n);
    copy_matrix(m, n, &A_tmp[0], m, &A[0], lda);

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[i + j*lda] = A_tmp[p[i] + j*m];
        }
    }
}

template <typename T>
void permute_rows_inv( const int seed, const int m, const int n, T *A, const int lda ){
    std::vector<int> p(m);
    for ( int i = 0; i < m; ++i ){
        p[i] = i;
    }
    std::mt19937 gen(seed);
    std::shuffle(p.begin(), p.end(), gen);

    std::vector<T> A_tmp(m * n);
    copy_matrix(m, n, &A_tmp[0], m, &A[0], lda);

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[p[i] + j*lda] = A_tmp[i + j*m];
        }
    }
}

template <typename T>
void permute_cols( const int seed, const int m, const int n, T *A, const int lda ){
    std::vector<int> p(n);
    for ( int i = 0; i < n; ++i ){
        p[i] = i;
    }
    std::mt19937 gen(seed);
    std::shuffle(p.begin(), p.end(), gen);

    std::vector<T> A_tmp(m * n);
    copy_matrix(m, n, &A_tmp[0], m, &A[0], lda);

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[i + j*lda] = A_tmp[i + p[j]*m];
        }
    }
}

template <typename T>
void permute_cols_inv( const int seed, const int m, const int n, T *A, const int lda ){
    std::vector<int> p(n);
    for ( int i = 0; i < n; ++i ){
        p[i] = i;
    }
    std::mt19937 gen(seed);
    std::shuffle(p.begin(), p.end(), gen);

    std::vector<T> A_tmp(m * n);
    copy_matrix(m, n, &A_tmp[0], m, &A[0], lda);

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[i + p[j]*lda] = A_tmp[i + j*m];
        }
    }
}

template <typename T>
void print_matrix_to_stdout_hex( const int m, const int n, const T *A, const int lda ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            printf("%a ", A[i+j*lda]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void print_matrix_to_stdout_sci( const int m, const int n, const T *A, const int lda ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            printf("%E ", A[i+j*lda]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void display_distribution_stats( const int m, const int n, const T *errors, const int lda, std::string name ) {

    std::vector<T> vals;
    vals.reserve(m * n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
            vals.push_back(errors[i + j * lda]);

    std::sort(vals.begin(), vals.end());
    const int count = vals.size();

    // Compute statistics
    auto percentile = [&](double p) -> T {
        double idx = p * (count - 1);
        int lo = (int)idx;
        int hi = lo + 1 < count ? lo + 1 : lo;
        double frac = idx - lo;
        return vals[lo] * (1.0 - frac) + vals[hi] * frac;
    };

    T p0   = vals.front();       // min
    T p25  = percentile(0.25);   // Q1
    T p50  = percentile(0.50);   // median
    T p75  = percentile(0.75);   // Q3
    T p100 = vals.back();        // max

    std::cout << "\n=== Distribution Stats | " << name << " ===\n\n";
    std::cout << "  min    : " << p0   << "\n";
    std::cout << "  Q1     : " << p25  << "\n";
    std::cout << "  median : " << p50  << "\n";
    std::cout << "  Q3     : " << p75  << "\n";
    std::cout << "  max    : " << p100 << "\n";
    std::cout << "  count  : " << count << "\n\n";
}

#endif