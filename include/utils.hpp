#ifndef UTILS_HPP
#define UTILS_HPP

#ifndef DRYRUN
#include <H5Cpp.h>

#define METADATA_STR_LEN 128

using namespace H5;

struct GEMM_Metadata {
    int matrix_id;
    char matrix_name[METADATA_STR_LEN];
    char input_type[METADATA_STR_LEN];
    double scaling_factor;
    int m;
    int n;
    int seed;
};
#endif

using namespace std;


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

#ifndef DRYRUN
H5File init_hdf5( string filename ){
    
    H5::Exception::dontPrint();
    H5File h5_file;
    try{
        h5_file = H5File(filename, H5F_ACC_RDWR);
    } catch ( H5::FileIException& e ){
        h5_file = H5File(filename, H5F_ACC_TRUNC);
        Group matrices_group;
        string matrices_group_name = "matrices";
        matrices_group = h5_file.createGroup(matrices_group_name);

        H5::CompType meta_type(sizeof(GEMM_Metadata));

        meta_type.insertMember(
            "matrix_id",
            HOFFSET(GEMM_Metadata, matrix_id),
            H5::PredType::NATIVE_INT32
        );

        meta_type.insertMember(
            "matrix_name",
            HOFFSET(GEMM_Metadata, matrix_name),
            H5::StrType(H5::PredType::C_S1, METADATA_STR_LEN)
        );

        meta_type.insertMember(
            "input_type",
            HOFFSET(GEMM_Metadata, input_type),
            H5::StrType(H5::PredType::C_S1, METADATA_STR_LEN)
        );

        meta_type.insertMember(
            "m",
            HOFFSET(GEMM_Metadata, m),
            H5::PredType::NATIVE_INT32
        );

        meta_type.insertMember(
            "n",
            HOFFSET(GEMM_Metadata, n),
            H5::PredType::NATIVE_INT32
        );

        meta_type.insertMember(
            "seed",
            HOFFSET(GEMM_Metadata, seed),
            H5::PredType::NATIVE_INT32
        );

        meta_type.insertMember(
            "scaling_factor",
            HOFFSET(GEMM_Metadata, scaling_factor),
            H5::PredType::NATIVE_DOUBLE
        );

        hsize_t initial_size = 0;
        hsize_t max_size     = H5S_UNLIMITED;

        H5::DataSpace meta_space(1, &initial_size, &max_size);

        // Enable chunking (required for extendable datasets)
        H5::DSetCreatPropList plist;
        hsize_t meta_chunk_size = 1024;
        plist.setChunk(1, &meta_chunk_size);
        plist.setDeflate(4);

        H5::DataSet meta_ds = h5_file.createDataSet(
            "/metadata",
            meta_type,
            meta_space,
            plist
        );
    }
    return h5_file;
}

int construct_matrix_id( string matrix_name, string input_type, double scaling_factor, int m, int n, int seed ){
    string id_str = matrix_name + "_" + input_type + "_" + to_string(scaling_factor) + "_" + to_string(m) + "m_" + to_string(n) + "n_" + to_string(seed);
    std::hash<string> hasher;
    return static_cast<int>( hasher(id_str) );
}

template <typename T> const PredType & native_type;
template <> const PredType& native_type<double> = PredType::NATIVE_DOUBLE;
template <> const PredType& native_type<long double> = PredType::NATIVE_LDOUBLE;

template <typename T>
void write_matrix_to_hdf5( H5File& h5_file, T *A, const int lda, string matrix_name, string input_type, double scaling_factor, int m, int n, int seed ){
    
    string matrices_group_name = "matrices";
    Group matrices_group = h5_file.openGroup(matrices_group_name);
    int matrix_id = construct_matrix_id(matrix_name, input_type, scaling_factor, m, n, seed);
    string matrix_id_name = to_string(matrix_id);
    
    DSetCreatPropList plist;
    hsize_t matrix_chunk_size[2] = {(hsize_t)m, 1};
    plist.setChunk(2, matrix_chunk_size);
    plist.setDeflate(4);
    
    hsize_t dimsA[2] = { (hsize_t)m, (hsize_t)n };
    DataSpace A_file_space(2, dimsA);
    DataSet dset;
    if ( H5Lexists(matrices_group.getId(), matrix_id_name.c_str(), H5P_DEFAULT) ){
        if ( (matrix_name == "A") || (matrix_name == "B") ){
            return;
        }
        dset = matrices_group.openDataSet(matrix_id_name.c_str(), H5P_DEFAULT);
    }
    else{
        dset = matrices_group.createDataSet(matrix_id_name, native_type<T>, A_file_space, plist);
        DataSet metadata_ds = h5_file.openDataSet("/metadata");
        DataSpace metadata_dataspace_old = metadata_ds.getSpace();
        hsize_t old_size;
        metadata_dataspace_old.getSimpleExtentDims(&old_size);

        hsize_t new_size = old_size + 1;
        metadata_ds.extend(&new_size);

        H5::DataSpace metadata_file_space = metadata_ds.getSpace();
        hsize_t offset = old_size;
        hsize_t count  = 1;
        metadata_file_space.selectHyperslab(H5S_SELECT_SET, &count, &offset);

        // Memory space
        H5::DataSpace metadata_mem_space(1, &count);

        // Fill metadata record
        GEMM_Metadata metadata{};
        metadata.matrix_id = matrix_id;
        snprintf(metadata.matrix_name, METADATA_STR_LEN, "%s", matrix_name.c_str());
        snprintf(metadata.input_type, METADATA_STR_LEN, "%s", input_type.c_str());
        metadata.scaling_factor = scaling_factor;
        metadata.m  = m;
        metadata.n  = n;
        metadata.seed = seed;

        // Write
        H5::CompType meta_type = metadata_ds.getCompType();
        metadata_ds.write(&metadata, meta_type, metadata_mem_space, metadata_file_space);
    }

    for ( int j = 0; j < n; ++j ){
        hsize_t count[2] = { (hsize_t)n, 1 };
        hsize_t start[2] = { 0, (hsize_t)j };
        A_file_space.selectHyperslab( H5S_SELECT_SET, count, start );

        hsize_t mem_dims[1] = { (hsize_t)n };
        DataSpace A_mem_space(1, mem_dims);

        dset.write(
            A + j * lda,
            native_type<T>,
            A_mem_space,
            A_file_space
        );
    }
}
#endif
#endif