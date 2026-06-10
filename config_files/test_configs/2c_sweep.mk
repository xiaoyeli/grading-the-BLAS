# === test parameters ===
INPUT_TYPE := test2c

LOG2_MIN_DIM := 5				# for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
LOG2_MAX_DIM := 10				#       // A and B have shape (2^d) x (2^d)
                                #       ...
LOG2_MIN_SCALING_FACTOR := 0    #       for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
LOG2_MAX_SCALING_FACTOR := 500  #           // entries in A and B span the orders of 2^-s to 2^s (note that LOG2_INC_SCALING_FACTOR must be >0)
LOG2_INC_SCALING_FACTOR := 100  #           ...
                                #           // execute DGEMM under test
                                #           ...
                                # }

# === advanced options ===
REFERENCE_GEMM := 2 # when REFERENCE_GEMM=1, reference is a naive triple loop GEMM in double precision
                    # when REFERENCE_GEMM=2, reference is a naive triple loop GEMM in long double precision

DUMP_A_B := 0       # when DUMP_A_B=1, dump the binary data of the input matrices to results/raw_data

DUMP_C := 0         # when DUMP_C=1, dump the binary data of the output matrices to results/raw_data

VERBOSE := 0        # when VERBOSE=1, prints A, B, and C to stdout

OUTER_PERMUTE := 0  # for OUTER_PERMUTE=n: Perform n trials of P1_inv * ( (P1 * A) * (B * P2) ) * P2_inv
                    # with random P1 and P2. The value of REFERENCE_GEMM is
                    # ignored and the maximum componentwise error with respect
                    # to A * B using the implementation under test over the n
                    # trials is reported

INNER_PERMUTE := 0  # for INNER_PERMUTE=n: Perform n trials of (A * P) * (P * B)
                    # with random P. The value of REFERENCE_GEMM is
                    # ignored and the maximum componentwise error with respect
                    # to A * B using the implementation under test over the n
                    # trials is reported