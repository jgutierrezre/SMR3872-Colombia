/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz = 0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms, nfi, nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    double *cx, *cy, *cz;
    int nsize, mpirank;
};
typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE* fp, char* buf) {
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp, BLEN, fp)) {
        int i;

        ptr = strchr(tmp, '#');
        if (ptr)
            *ptr = '\0';
        i = strlen(tmp);
        --i;
        while (isspace(tmp[i])) {
            tmp[i] = '\0';
            --i;
        }
        ptr = tmp;
        while (isspace(*ptr)) {
            ++ptr;
        }
        i = strlen(ptr);
        strcpy(buf, tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}

/* helper function: get current time in seconds since epoch */

static double wallclock() {
    struct timeval t;
    gettimeofday(&t, 0);
    return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

/* helper function: zero out an array */
static void azzero(double* d, const int n) {
    int i;
    for (i = 0; i < n; ++i) {
        d[i] = 0.0;
    }
}

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2) {
    while (x > boxby2)
        x -= 2.0 * boxby2;
    while (x < -boxby2)
        x += 2.0 * boxby2;
    return x;
}

/* compute kinetic energy */
static void ekin(mdsys_t* sys) {
    int i;

    sys->ekin = 0.0;
    for (i = 0; i < sys->natoms; ++i) {
        sys->ekin += 0.5 * mvsq2e * sys->mass *
                     (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] +
                      sys->vz[i] * sys->vz[i]);
    }
    sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}

/* compute forces */
static void force(mdsys_t* sys) {
    double rsq, ffac;
    double rx, ry, rz;
    double c6, c12, rcsq;
    int i, j, ii;

    /* zero energy and forces */
    sys->epot = 0.0;
    double epot = 0.0;
    azzero(sys->cx, sys->natoms);
    azzero(sys->cy, sys->natoms);
    azzero(sys->cz, sys->natoms);
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Datatype mpi_sys_type;
    // MPI_Bcast( sys , 1 , mpi_sys_type, 0 , MPI_COMM_WORLD);
    // MPI_Comm_rank(MPI_COMM_WORLD, &(sys->mpirank));

    c6 = pow(sys->sigma, 6.0);
    c12 = 4.0 * sys->epsilon * c6 * c6;
    c6 *= 4.0 * sys->epsilon;
    rcsq = sys->rcut * sys->rcut;
    for (i = 0; i < (sys->natoms - 1);
         i += sys->nsize) {  // In the slides it goes up to sys->natoms-1, but I
                             // think it should be sys->natoms
        ii = i + sys->mpirank;
        if (ii >= (sys->natoms - 1))
            break;
        for (j = ii + 1; j < (sys->natoms);
             ++j) {  // I'm using Newton's 3rd law here

            /* get distance between particle i and j */
            rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
            rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
            rsq = rx * rx + ry * ry + rz * rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double r6, rinv;
                rinv = 1.0 / rsq;
                r6 = rinv * rinv * rinv;
                ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;

                epot += r6 * (c12 * r6 - c6);

                sys->cx[ii] += rx * ffac;
                sys->cy[ii] += ry * ffac;
                sys->cz[ii] += rz * ffac;
                // Newton's 3rd law
                sys->cx[j] -= rx * ffac;
                sys->cy[j] -= ry * ffac;
                sys->cz[j] -= rz * ffac;
            }
        }
    }
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/* velocity verlet */
static void velverlet(mdsys_t* sys) {
    int i;

    if (sys->mpirank == 0) {
        /* first part: propagate velocities by half and positions by full step
         */
        for (i = 0; i < sys->natoms; ++i) {
            sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
            sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
            sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
            sys->rx[i] += sys->dt * sys->vx[i];
            sys->ry[i] += sys->dt * sys->vy[i];
            sys->rz[i] += sys->dt * sys->vz[i];
        }
    }
    /*    MPI_Barrier( MPI_COMM_WORLD);
        MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->vx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->vy, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->vz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
    // printf("I'm %d and I'm here\n", sys->mpirank);
    //    MPI_Barrier( MPI_COMM_WORLD);
    /* compute forces and potential energy */
    force(sys);
    /* second part: propagate velocities by another half step */
    if (sys->mpirank == 0) {
        for (i = 0; i < sys->natoms; ++i) {
            sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
            sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
            sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        }
    }
    /*    MPI_Barrier( MPI_COMM_WORLD);
        MPI_Bcast(sys->vx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->vy, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->vz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier( MPI_COMM_WORLD);*/
}

/* append data to output. */
static void output(mdsys_t* sys, FILE* erg, FILE* traj) {
    int i;

    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
           sys->ekin, sys->epot, sys->ekin + sys->epot);
    fprintf(erg, "% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
            sys->ekin, sys->epot, sys->ekin + sys->epot);
    fprintf(traj, "%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi,
            sys->ekin + sys->epot);
    for (i = 0; i < sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i],
                sys->rz[i]);
    }
}

/* main */
int main(int argc, char** argv) {
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp, *traj, *erg;
    mdsys_t sys;
    double t_start;
    /* MPI initialization */
    MPI_Init(&argc, &argv);
    // MPI_COMM_WORLD=MPI_COMM_WORLD;
    MPI_Comm_size(MPI_COMM_WORLD, &sys.nsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &sys.mpirank);

    if (sys.mpirank == 0) {
        printf("LJMD version %3.1f\n", LJMD_VERSION);

        t_start = wallclock();

        /* read input file */
        if (get_a_line(stdin, line))
            return 1;
        sys.natoms = atoi(line);
        if (get_a_line(stdin, line))
            return 1;
        sys.mass = atof(line);
        if (get_a_line(stdin, line))
            return 1;
        sys.epsilon = atof(line);
        if (get_a_line(stdin, line))
            return 1;
        sys.sigma = atof(line);
        if (get_a_line(stdin, line))
            return 1;
        sys.rcut = atof(line);
        if (get_a_line(stdin, line))
            return 1;
        sys.box = atof(line);
        if (get_a_line(stdin, restfile))
            return 1;
        if (get_a_line(stdin, trajfile))
            return 1;
        if (get_a_line(stdin, ergfile))
            return 1;
        if (get_a_line(stdin, line))
            return 1;
        sys.nsteps = atoi(line);
        if (get_a_line(stdin, line))
            return 1;
        sys.dt = atof(line);
        if (get_a_line(stdin, line))
            return 1;
        nprint = atoi(line);
    }
    // printf("natoms=%d\n", sys.natoms);
    int block_lengths[18] = {1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[18] = {
        MPI_INT,    MPI_INT,    MPI_INT,    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Aint offsets[18];

    // Calculate offsets for each member
    offsets[0] = offsetof(mdsys_t, natoms);
    offsets[1] = offsetof(mdsys_t, nfi);
    offsets[2] = offsetof(mdsys_t, nsteps);
    offsets[3] = offsetof(mdsys_t, dt);
    offsets[4] = offsetof(mdsys_t, mass);
    offsets[5] = offsetof(mdsys_t, epsilon);
    offsets[6] = offsetof(mdsys_t, sigma);
    offsets[7] = offsetof(mdsys_t, box);
    offsets[8] = offsetof(mdsys_t, rcut);
    offsets[9] = offsetof(mdsys_t, ekin);
    offsets[10] = offsetof(mdsys_t, epot);
    offsets[11] = offsetof(mdsys_t, temp);
    offsets[12] = offsetof(mdsys_t, rx);
    offsets[13] = offsetof(mdsys_t, ry);
    offsets[14] = offsetof(mdsys_t, rz);
    offsets[15] = offsetof(mdsys_t, vx);
    offsets[16] = offsetof(mdsys_t, vy);
    offsets[17] = offsetof(mdsys_t, vz);

    // Create the MPI datatype
    MPI_Datatype mpi_sys_type;
    MPI_Type_create_struct(18, block_lengths, offsets, types, &mpi_sys_type);
    MPI_Type_commit(&mpi_sys_type);

    // MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&nprint, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys, 1, mpi_sys_type, 0, MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &sys.mpirank);
    // printf("natoms=%d\n", sys.natoms);
    /* allocate memory */
    sys.rx = (double*)malloc(sys.natoms * sizeof(double));
    sys.ry = (double*)malloc(sys.natoms * sizeof(double));
    sys.rz = (double*)malloc(sys.natoms * sizeof(double));
    sys.vx = (double*)malloc(sys.natoms * sizeof(double));
    sys.vy = (double*)malloc(sys.natoms * sizeof(double));
    sys.vz = (double*)malloc(sys.natoms * sizeof(double));
    sys.fx = (double*)malloc(sys.natoms * sizeof(double));
    sys.fy = (double*)malloc(sys.natoms * sizeof(double));
    sys.fz = (double*)malloc(sys.natoms * sizeof(double));
    sys.cx = (double*)malloc(sys.natoms * sizeof(double));
    sys.cy = (double*)malloc(sys.natoms * sizeof(double));
    sys.cz = (double*)malloc(sys.natoms * sizeof(double));

    if (sys.mpirank == 0) {
        /* read restart */
        fp = fopen(restfile, "r");
        if (fp) {
            for (i = 0; i < sys.natoms; ++i) {
                fscanf(fp, "%lf%lf%lf", sys.rx + i, sys.ry + i, sys.rz + i);
            }
            for (i = 0; i < sys.natoms; ++i) {
                fscanf(fp, "%lf%lf%lf", sys.vx + i, sys.vy + i, sys.vz + i);
            }
            fclose(fp);
            azzero(sys.fx, sys.natoms);
            azzero(sys.fy, sys.natoms);
            azzero(sys.fz, sys.natoms);
        } else {
            perror("cannot read restart file");
            return 3;
        }
    }
    /* initialize forces and energies.*/
    //    MPI_Barrier( MPI_COMM_WORLD);
    // printf("Initialize forces and energies\n");
    sys.nfi = 0;

    force(&sys);
    // printf("Computed first force\n");
    if (sys.mpirank == 0) {
        ekin(&sys);
        erg = fopen(ergfile, "w");
        traj = fopen(trajfile, "w");

        printf("Startup time: %10.3fs\n", wallclock() - t_start);
        printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms,
               sys.nsteps);
        printf(
            "     NFI            TEMP            EKIN                 EPOT     "
            "         ETOT\n");
        output(&sys, erg, traj);

        /* reset timer */
        t_start = wallclock();
    }
    /**************************************************/
    /* main MD loop */
    for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0 && sys.mpirank == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet(&sys);
        if (sys.mpirank == 0)
            ekin(&sys);
    }
    /**************************************************/

    if (sys.mpirank == 0) {
        /* clean up: close files, free memory */
        printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
        fclose(erg);
        fclose(traj);
    }

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    free(sys.cx);
    free(sys.cy);
    free(sys.cz);

    MPI_Finalize();
    return 0;
}