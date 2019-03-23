#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"
#include<iostream>
#include<vector>

#define NUM_THREADS 256



struct Bin{
    particle_t** particles;
    int number_of_particles;
};




extern double size;
//
//  benchmarking program
//

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}


//original code
/* 
__global__ void compute_forces_gpu(particle_t * particles, int n)
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  for(int j = 0 ; j < n ; j++)
    apply_force_gpu(particles[tid], particles[j]);

}
*/



__global__ void compute_forces_gpu(particle_t * particles, int n)
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  for(int j = 0 ; j < n ; j++)
    apply_force_gpu(particles[tid], particles[j]);

}






__global__ void move_gpu (particle_t * particles, int n, double size)
{

  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}


//method called by host to send the grid to gpu
__host__ Bin* send_grid_to_gpu(Bin* grid, int dim){
    Bin* d_grid;
    cudaMalloc((void **) &d_grid, dim * dim * sizeof(Bin));
    cudaMemcpy(d_grid, grid, dim * dim * sizeof(Bin), cudaMemcpyHostToDevice);
    return d_grid;
}


__device__ void apply_forces_to_cell(Bin &src, Bin &neighbor){
    int num_particles_src = src.number_of_particles;
    int num_particles_neighbor = neighbor.number_of_particles;

    for(int i = 0; i < num_particles_src; i++){
        for(int j = 0; j < num_particles_neighbor; j++){
            apply_force_gpu(*(src.particles[i]), *(neighbor.particles[j]));
        }
    }

}

//method to apply the forces on the bins
__device__ void apply_forces_on_grid(Bin* grid, const int dim, int tid){

    // assume that # of threads are <dim>
    int bin_i = tid/dim;
    int bin_j = tid%dim;
    int index = bin_i * dim + bin_j;
    int number_of_particles = grid[index].number_of_particles;
    
    // initilize acceleration
    for(int i = 0; i < number_of_particles; i++){
        grid[index].particles[i]->ax = 0;
        grid[index].particles[i]->ay = 0;        
    }
    
    for(int i = -1; i < 2; i++){
        int delta_i = bin_i + i;
        for(int j = -1; j < 2; j++){
            int delta_j = bin_j + j;
            if(delta_i >= 0 && delta_j >= 0 && delta_i < dim && delta_j < dim){
                int index2 = delta_i * dim + delta_j;
                apply_forces_to_cell(grid[index], grid[index2]);
            }
        }
    }
   
} 



int main( int argc, char **argv )
{    
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    set_size( n );

    init_particles( n, particles );


    std::cout << "size: " << size << std::endl;
    double world_dim = size;
    int grid_dim = int(ceil(sqrt(n)));
    double cell_size = world_dim/grid_dim;
    //ensure to not violate cutoff constraint
    // if violated, set cell_size to be the cutoff
    std::cout<< "cutoff: " << cutoff << std::endl; 
    if(cell_size < cutoff){
        grid_dim = int(world_dim / cutoff);
        cell_size= world_dim / grid_dim; 
    }

    //init grid
    Bin* grid = (Bin*) malloc(grid_dim*grid_dim*sizeof(Bin));
    for(int i = 0; i < grid_dim*grid_dim; i++)
        grid[i] = Bin(); 

    //used vector in host temporarily to bin particles in host
    // then papck the particles in the grid array to be used in gpu
    std::vector<particle_t*> temp[grid_dim][grid_dim];
    for(int i = 0; i < n; i++){
        int index_i = floor(particles[i].y/cell_size);
        int index_j = floor(particles[i].x/cell_size);
        temp[index_i][index_j].push_back(&particles[i]);
    }
    //naive way,, improve later
    for(int i = 0; i < grid_dim; i++){
        for(int j = 0; j < grid_dim; j++){
            int number_of_particles = temp[i][j].size();
            grid[i*grid_dim + j].number_of_particles = number_of_particles;
            grid[i*grid_dim + j].particles = (particle_t**) malloc(number_of_particles * sizeof(particle_t*));
            for(int k = 0; k < number_of_particles; k++){
                grid[i*grid_dim + j].particles[k] = temp[i][j][k];
            }
            temp[i][j].clear();     
        }
    }

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //

	int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
	compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n);
        
        //
        //  move particles
        //
	move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
	    // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
	}
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
