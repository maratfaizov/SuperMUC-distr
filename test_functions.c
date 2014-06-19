/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util_read_files.h"
#include "util_write_files.h"
#include "mpi.h"


int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *cgup) {
    int NINTCI, NINTCF, NEXTCI, NEXTCF, points_count;
    int **LCC, **points;
    double *BS, *BE, *BN, *BW, *BL, *BH;
    double *BP, *SU, *distr;
    int *elems;

    int i, my_rank, num_procs, i_end;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int CHUNKSIZE =  (int)ceil((double)local_num_elems / (double)num_procs);  // ceil
    int REMAINING = (local_num_elems) - (num_procs-1) * CHUNKSIZE;

    printf("local_num_elems: %d \n", local_num_elems);

    i_end = CHUNKSIZE;

    if (my_rank == num_procs - 1) {
        i_end = REMAINING;
    }

    int f_status = read_binary_geo_single( file_in, &NINTCI, &NINTCF, &NEXTCI, &NEXTCF,
                                     &LCC, &BS, &BE, &BN, &BW, &BL, &BH,
                                           &BP, &SU, &points_count, &points, &elems);

    if ( f_status != 0 ) return f_status;

printf("nintcf: %d \n", NINTCF);
printf("nextcf: %d \n", NEXTCF);

    distr = (double*) calloc(sizeof(double), local_num_elems);
    printf("local_global_index[100]=%d \n", local_global_index[100]);
    printf("cgup[100]=%d \n", cgup[100]);
    printf("cgup[i_end-1]=%d \n", cgup[i_end-1]);
    printf("local_global_index[103]=%d \n", local_global_index[i_end-1]);
    for ( i = 0 ; i < i_end; i++ ) {
        distr[local_global_index[i]]=cgup[i];
    }

    printf("TEST after read 2 \n");

    vtk_write_unstr_grid_header( "distribution_test", file_vtk_out,
                                 NINTCI, NINTCF, points_count, points, elems);

    vtk_append_double(file_vtk_out, "CGUP", 0, NINTCF, distr);

    free(SU);
    free(BP);
    free(BH);
    free(BL);
    free(BW);
    free(BN);
    free(BE);
    free(BS);
    free(elems);

    printf("test_function succesfully terminated\n");

    return 0;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index,
                       int local_num_elems, int neighbors_count, int* send_count, int** send_list,
                       int* recv_count, int** recv_list) {
    // Return an error if not implemented
    return -1;
}

