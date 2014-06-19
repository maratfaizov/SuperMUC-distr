/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util_read_files.h"
#include "initialization.h"
#include "mpi.h"
#include "metis.h"


int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
        double** bl, double** bh, double** bp, double** su, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
        int*** recv_list, int** epart, int** npart, int* objval) {
    /********** START INITIALIZATION **********/


    int f_status = 0;

    /********** METIS DISTRIBUTION ************/

    if ( ( strcmp( part_type, "dual" ) == 0 ) || ( strcmp( part_type, "nodal" ) == 0 ) ) {
        f_status = initialization_metis(file_in, part_type, nintci, nintcf, nextci,
                nextcf, lcc, bs, be, bn, bw,
                bl, bh, bp, su, points_count,
                points, elems, var, cgup, oc,
                cnorm, local_global_index, global_local_index,
                neighbors_count, send_count, send_list, recv_count,
                recv_list, epart, npart, objval);
    } else {
    /******** CLASSIC DISTRIBUTION ************/
        f_status = initialization_classic(file_in, part_type, nintci, nintcf, nextci,
                nextcf, lcc, bs, be, bn, bw,
                bl, bh, bp, su, points_count,
                points, elems, var, cgup, oc,
                cnorm, local_global_index, global_local_index,
                neighbors_count, send_count, send_list, recv_count,
                recv_list, epart, npart, objval);
    }

    return f_status;
}



int initialization_classic(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
        double** bl, double** bh, double** bp, double** su, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
        int*** recv_list, int** epart, int** npart, int* objval) {
    int i = 0;
    int my_rank, num_procs;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    // read-in the input file

    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
            &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
            &*points, &*elems);

    if ( f_status != 0 ) return f_status;

    // size of internal and external junks
    int CHUNKSIZE_INT =  (int)ceil((double)(*nintcf-*nintci+1) / (double)num_procs);  // ceil
    int REMAINING_INT = (*nintcf-*nintci +1) - (num_procs-1) * CHUNKSIZE_INT;

    int CHUNKSIZE_EXT =  (int)ceil((double)(*nextcf-*nextci+1) / (double)num_procs);  // ceil
    int REMAINING_EXT = (*nextcf-*nextci +1) - (num_procs-1) * CHUNKSIZE_EXT;

    *var = (double*) calloc(sizeof(double), (CHUNKSIZE_INT + CHUNKSIZE_EXT));
    *cgup = (double*) calloc(sizeof(double), (CHUNKSIZE_INT + CHUNKSIZE_EXT));
    *cnorm = (double*) calloc(sizeof(double), (CHUNKSIZE_INT));

    int i_end_int = CHUNKSIZE_INT;
    int i_end_ext = CHUNKSIZE_EXT;

    if (my_rank == num_procs - 1) {
        i_end_int = REMAINING_INT;
        i_end_ext = REMAINING_EXT;
    }

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        (*cnorm)[i] = 1.0;
    }

    for ( i = 0; i < i_end_int; i++ ) {
        (*var)[i] = 0.0;
    }

    for ( i = 0; i < i_end_int; i++ ) {
        (*cgup)[i] = 1.0 / ((*bp)[i]);
    }
    printf("i_end_int=%d , i_end_ext=%d, CHUNKSIZE_INT = %d", i_end_int, i_end_ext, CHUNKSIZE_INT);
MPI_Barrier( MPI_COMM_WORLD );
    for ( i = CHUNKSIZE_INT; i < CHUNKSIZE_INT + i_end_ext; i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bh)[i] = 0.0;
        (*bl)[i] = 0.0;
    }

    // initialize local to global mapping

    *local_global_index = (int*) malloc(sizeof(int) * (CHUNKSIZE_INT+CHUNKSIZE_EXT));

    for ( i = 0; i < CHUNKSIZE_INT; i++ ) {
        (*local_global_index)[i] = i + my_rank * CHUNKSIZE_INT;
    }

    for ( i = CHUNKSIZE_INT; i < CHUNKSIZE_INT + CHUNKSIZE_EXT; i++ ) {
        (*local_global_index)[i] = *nintcf-*nintci + i + my_rank * CHUNKSIZE_EXT;
    }

    return 0;
}



int initialization_metis(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
        double** bl, double** bh, double** bp, double** su, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
        int*** recv_list, int** epart, int** npart, int* objval) {
    int i = 0, j = 0;
    int my_rank, num_procs, CHUNKSIZE, CHUNKSIZE_TOT;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // define CHUNKSIZE

    printf("nodal or dual\n");
    int CHUNKSIZE0 = (int)ceil((double)(*nintcf - *nintci + 1 ) / (double) num_procs);
    int REMAINDER = ( *nintcf - *nintci + 1 ) - ( num_procs - 1 ) * CHUNKSIZE0;

    if (my_rank == num_procs -1) {
        CHUNKSIZE = REMAINDER;
    } else {
        CHUNKSIZE = CHUNKSIZE0;
    }

    int CHUNKSIZE0_TOT = (int)ceil((double)(*nextcf - *nintci + 1 ) / (double) num_procs);
    int REMAINDER_TOT = ( *nextcf - *nintci + 1 ) - ( num_procs - 1 ) * CHUNKSIZE0_TOT;

    if (my_rank == num_procs -1) {
        CHUNKSIZE_TOT = REMAINDER_TOT;
    } else {
        CHUNKSIZE_TOT = CHUNKSIZE0_TOT;
    }
    printf("TEST1\n");
    MPI_Barrier( MPI_COMM_WORLD);

    // initialize global data structures

    double *gbs, *gbe, *gbn, *gbw, *gbl, *gbh, *gbp, *gsu;
    int gpoints_count;
    int** gpoints;
    int* gelems;

    double *gvar, *gcgup, *gcnorm;


    //    printf("Process %i initializated function/n", my_rank);

    // read-in the input file
    int f_status = read_binary_geo_single( file_in, nintci, nintcf, nextci,
            nextcf, lcc, &gbs, &gbe, &gbn,
            &gbw, &gbl, &gbh, &gbp, &gsu,
            &gpoints_count,    &gpoints, &gelems );

    printf( "rank %d: binary read in successful code %d \n", my_rank, f_status );
    MPI_Barrier( MPI_COMM_WORLD );

    if ( f_status != 0 ) return f_status;

    gvar = (double*) calloc( (*nextcf + 1), sizeof(double) );
    gcgup = (double*) calloc( (*nextcf + 1), sizeof(double) );
    gcnorm = (double*) calloc( (*nintcf + 1), sizeof(double) );

    // initialise the arrays
    for ( i = 0; i <= 10; i++ ) {
        gcnorm[i] = 1.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        gvar[i] = 0.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        gcgup[i] = 1.0 / ((gbp)[i]);
    }

    for ( i = (*nextci); i <= (*nextcf); i++ ) {
        gvar[i] = 0.0;
        gcgup[i] = 0.0;
        gbs[i] = 0.0;
        gbe[i] = 0.0;
        gbn[i] = 0.0;
        gbw[i] = 0.0;
        gbh[i] = 0.0;
        gbl[i] = 0.0;
    }

    /************METIS************/

    idx_t ne;
    idx_t nn;
    idx_t *eptr;
    idx_t *eind;
    idx_t ncommon;
    idx_t nparts;
    idx_t *idx_objval;
    idx_t *idx_epart;
    idx_t *idx_npart;

    ne = ( *nintcf )  - ( *nintci ) + 1;
    nn = gpoints_count;

    ncommon = 4;
    nparts = num_procs;
    printf("TEST_1\n");
    MPI_Barrier( MPI_COMM_WORLD);

    eptr = (idx_t*) calloc( ( ne + 1 ), sizeof(idx_t) );
    eind = (idx_t*) calloc( ( ( ne + 1 ) * 8 ), sizeof(idx_t) );

    idx_objval = (idx_t*) calloc( 1, sizeof(idx_t) );
    idx_epart = (idx_t*) calloc( ( ne ), sizeof(idx_t) );
    idx_npart = (idx_t*) calloc( ( nn ), sizeof(idx_t) );

    epart = (int**) calloc( ( ne ) , sizeof(int*) );
    npart = (int**) calloc( ( nn ),  sizeof(int*) );

    for ( i = 0; i < ( ne + 1 ); i++ ) {
        eptr[i] = i * 8;
    }

    for ( i = 0; i < ( ( ne + 1 ) * 8 ); i++ ) {
        eind[i] = (gelems)[i];
    }

    if ( strcmp( part_type, "dual" ) == 0 ) {
        METIS_PartMeshDual( &ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts,
                NULL, NULL, idx_objval, idx_epart, idx_npart );
    } else {
        METIS_PartMeshNodal( &ne, &nn,    eptr, eind, NULL, NULL, &nparts,
                NULL, NULL, idx_objval, idx_epart, idx_npart );
    }


    printf("idx_epart[1000]=%d \n", idx_epart[ne-1]);
    MPI_Barrier( MPI_COMM_WORLD);


    for (i = 0; i < ne; i++) {
        epart[i]=(int)idx_epart[i];
    }

    for (i = 0; i < nn; i++) {
        npart[i]=(int)idx_npart[i];
    }

    *objval=(int)*idx_objval;


    // local_global_index
    if ( ( *local_global_index = (int*) calloc( ne , sizeof(int) ) ) == NULL ) {
        fprintf(stderr, "malloc failed to allocate local_global_index\n");
        return -1;
    }

    if ( (  *global_local_index = (int*) calloc(ne, sizeof(int) ) ) == NULL ) {
        fprintf(stderr, "malloc failed to allocate local_global_index\n");
        return -1;
    }

    j = 0;
    for (i = 0; i < ne; i++) {
        if( epart[i] == my_rank ) {
            (*local_global_index)[j]= i;
                        (*global_local_index)[i]= j;
            j++;
        }
    }
        printf("local_global_index[100]=%d \n", local_global_index[100]);
        printf("global_local_index[100]=%d \n", global_local_index[100]);
        printf("epart[ne-1]=%d \n", epart[ne-1]);

        MPI_Barrier( MPI_COMM_WORLD);


    // allocate other arrays

    if ( (*cgup = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }

    if ( (*bs = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BS) failed\n");
        return -1;
    }

    if ( (*be = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BE) failed\n");
        return -1;
    }

    if ( (*bn = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BN) failed\n");
        return -1;
    }

    if ( (*bw = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BW) failed\n");
        return -1;
    }

    if ( (*bl = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BL) failed\n");
        return -1;
    }

    if ( (*bh = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BH) failed\n");
        return -1;
    }

    if ( (*bp = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BP) failed\n");
        return -1;
    }

    if ( (*su = (double *) malloc(CHUNKSIZE_TOT * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }


    printf("TEST_2\n");
    MPI_Barrier( MPI_COMM_WORLD);
    j = 0;
    for ( i = 0; i < ne; i++ ) {
        if(epart[i] == my_rank) {
            (*cgup)[j] = gcgup[i];

            (*bs)[j] = gbs[i];
            (*be)[j] = gbe[i];
            (*bn)[j] = gbn[i];
            (*bw)[j] = gbw[i];
            (*bl)[j] = gbl[i];
            (*bh)[j] = gbh[i];
            (*bp)[j] = gbp[i];
            (*su)[j] = gsu[i];
            j++;
        }
    }

    free(gbs);
    free(gbe);
    free(gbn);
    free(gbw);
    free(gbl);
    free(gbh);
    free(gbp);
    free(gsu);

    printf("TEST_3\n");
    MPI_Barrier( MPI_COMM_WORLD);


    return 0;
}
