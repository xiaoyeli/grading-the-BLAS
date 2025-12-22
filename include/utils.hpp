#ifndef UTILS_HPP
#define UTILS_HPP

#include <stdint.h>
#include <H5Cpp.h>

#define CHUNK_SIZE 32

using namespace H5;
using namespace std;

template <typename T>
void print_matrix_to_stdout_hex( const int32_t m, const int32_t n, const T *A, const int32_t lda ){
    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%a ", A[i+j*lda]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void print_matrix_to_stdout_dec( const int32_t m, const int32_t n, const T *A, const int32_t lda ){
    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%f ", A[i+j*lda]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}


void write_A_and_B_to_hdf5( H5File& h5_file, string matrix_type, int32_t seed, int32_t m, int32_t n, int32_t k, double *A, const int32_t lda, double *B, const int32_t ldb ){
    
    Group matrix_type_group;
    string matrix_type_group_name = "mat_" + matrix_type;
    if (H5Lexists(h5_file.getId(), matrix_type_group_name.c_str(), H5P_DEFAULT) ){
        matrix_type_group = h5_file.openGroup(matrix_type_group_name);
    } else {
        matrix_type_group = h5_file.createGroup(matrix_type_group_name);
    }
    
    Group dimension_group;
    string dimension_group_name = "dim_" + to_string(m) + "m_" + to_string(n) + "n_" + to_string(k) + "k";
    if ( H5Lexists(matrix_type_group.getId(), dimension_group_name.c_str(), H5P_DEFAULT) ){
        dimension_group = matrix_type_group.openGroup(dimension_group_name);
    } else {
        dimension_group = matrix_type_group.createGroup(dimension_group_name);
    }

    Group seed_group;
    string seed_group_name = "seed_" + to_string(seed);
    if ( H5Lexists(dimension_group.getId(), seed_group_name.c_str(), H5P_DEFAULT) ){
        seed_group = dimension_group.openGroup(seed_group_name);
    } else {
        seed_group = dimension_group.createGroup(seed_group_name);
    }


    if ( H5Lexists(seed_group.getId(), "A", H5P_DEFAULT) && H5Lexists(seed_group.getId(), "B", H5P_DEFAULT) ){
        return;
    }

    DSetCreatPropList plist;
    hsize_t chunk[2] = {CHUNK_SIZE, CHUNK_SIZE};
    plist.setChunk(2, chunk);
    plist.setDeflate(4);

    hsize_t dimsA[2] = { (hsize_t)m, (hsize_t)k };
    DataSpace file_spaceA(2, dimsA);
    DataSet dset = seed_group.createDataSet("A", PredType::NATIVE_DOUBLE, file_spaceA, plist);
    for ( int32_t j = 0; j < k; ++j ){
        hsize_t count[2] = { (hsize_t)k, 1 };
        hsize_t start[2] = { 0, (hsize_t)j };
        file_spaceA.selectHyperslab( H5S_SELECT_SET, count, start );

        hsize_t mem_dims[1] = { (hsize_t)k };
        DataSpace mem_spaceA(1, mem_dims);

        dset.write(
            A + j * lda,
            PredType::NATIVE_DOUBLE,
            mem_spaceA,
            file_spaceA
        );
    }

    hsize_t dimsB[2] = { (hsize_t)k, (hsize_t)n };
    DataSpace file_spaceB(2, dimsB);
    dset = seed_group.createDataSet("B", PredType::NATIVE_DOUBLE, file_spaceB, plist);
    for ( int32_t j = 0; j < n; ++j ){
        hsize_t start[2] = { 0, (hsize_t)j };
        hsize_t count[2] = { (hsize_t)n, 1 };
        file_spaceB.selectHyperslab( H5S_SELECT_SET, count, start );

        hsize_t mem_dims[1] = { (hsize_t)n };
        DataSpace mem_spaceB(1, mem_dims);

        dset.write(
            B + j * ldb,
            PredType::NATIVE_DOUBLE,
            mem_spaceB,
            file_spaceB
        );
    }
}

void write_C_to_hdf5( H5File& h5_file, string impl, string matrix_type, int32_t seed, int32_t m, int32_t n, int32_t k, double *C, const int32_t ldc){

    Group matrix_type_group;
    string matrix_type_group_name = "mat_" + matrix_type;
    if (H5Lexists(h5_file.getId(), matrix_type_group_name.c_str(), H5P_DEFAULT) ){
        matrix_type_group = h5_file.openGroup(matrix_type_group_name);
    } else {
        matrix_type_group = h5_file.createGroup(matrix_type_group_name);
    }
    
    Group dimension_group;
    string dimension_group_name = "dim_" + to_string(m) + "m_" + to_string(n) + "n_" + to_string(k) + "k";
    if ( H5Lexists(matrix_type_group.getId(), dimension_group_name.c_str(), H5P_DEFAULT) ){
        dimension_group = matrix_type_group.openGroup(dimension_group_name);
    } else {
        dimension_group = matrix_type_group.createGroup(dimension_group_name);
    }

    Group seed_group;
    string seed_group_name = "seed_" + to_string(seed);
    if ( H5Lexists(dimension_group.getId(), seed_group_name.c_str(), H5P_DEFAULT) ){
        seed_group = dimension_group.openGroup(seed_group_name);
    } else {
        seed_group = dimension_group.createGroup(seed_group_name);
    }

    Group impl_group;
    string impl_group_name = impl;
    if ( H5Lexists(seed_group.getId(), impl_group_name.c_str(), H5P_DEFAULT) ){
        impl_group = seed_group.openGroup(impl_group_name);
    } else{
        impl_group = seed_group.createGroup(impl_group_name);
    }

    hsize_t dims[2] = { (hsize_t)m, (hsize_t)n };
    DataSpace file_space(2, dims);
    DSetCreatPropList plist;
    hsize_t chunk[2] = {CHUNK_SIZE, CHUNK_SIZE};
    plist.setChunk(2, chunk);
    plist.setDeflate(4);
    DataSet dset = impl_group.createDataSet("C", PredType::NATIVE_DOUBLE, file_space, plist);

    for ( int32_t j = 0; j < n; ++j ){
        hsize_t count[2] = { (hsize_t)n, 1 };
        hsize_t start[2] = { 0, (hsize_t)j };
        file_space.selectHyperslab( H5S_SELECT_SET, count, start );

        hsize_t mem_dims[1] = { (hsize_t)n };
        DataSpace mem_space(1, mem_dims);

        dset.write(
            C + j * ldc,
            PredType::NATIVE_DOUBLE,
            mem_space,
            file_space
        );
    }

    return;
}
#endif