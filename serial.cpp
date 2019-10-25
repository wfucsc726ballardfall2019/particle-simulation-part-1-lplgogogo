#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <iostream>
using namespace std;

// constant copied from common.cpp for use
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //calculate the gridSize, binSize, and then number of bin on one side;
    double gridSize = sqrt(n * density);
    double binSize = cutoff * 2;     // equals to the diameter of the circle
    int binNum = int(gridSize / binSize) + 1; // the binNum should be +1
    cout << binNum <<endl;

    vector<vector<int> >bin(binNum * binNum, vector<int> (0));

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;

        //put all the particles into corresponding bins
        for (int i = 0; i < n; i++){
            int row = floor(particles[i].x / binSize);     //calculate the row index of the bin
            int col = floor(particles[i].y / binSize);     //calculate the column index of the bin
            bin[row * binNum + col].push_back(i);      //put the particle in to the bin in row major
        }

        for (int i = 0; i < n; i++){
            particles[i].ax = particles[i].ay = 0;      // initialize acceleration
            int row = floor(particles[i].x / binSize);     //calculate the row index of the bin
            int col = floor(particles[i].y / binSize);     //calculate the column index of the bin

            //
            //different cases to deal with
            //
            //situation that the particle is not in the first or the last column of the grid
            if ((row > 0) && (row< binNum-1)){
                for (int j = row-1; j <= row+1; j++){
                    //situation that the particle is not in the first or the last row of the grid
                    if ((col > 0 )&&(col < binNum-1)){
                        for (int k = col-1; k <= col+1; k++){
                             for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);
                            }
                        }
                    }

                    //situation that the particle is in the first column of the grid
                    else if (col == 0){
                        for (int k = col; k <= col+1; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);
                            }       
                        }
                    }
                    //situation that the particle is in the last column of the grid
                    else{
                        for (int k = col-1; k <= col; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                 
                            }
                        }
                    }
                }
            }

            //situation that the particle is in the first row of the grid
            else if (row == 0){
                for (int j = row; j <= row+1; j++){
                    //situation that the particle is not in the first or the last column of the grid
                    if ((col > 0 )&& (col < binNum-1)){
                        for (int k = col-1; k <= col+1; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                  
                            }
                        }
                    }
                    //situation that the particle is in the first column of the grid
                    else if (col == 0){
                        for (int k = col; k <= col+1; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                 
                            }
                        }
                    }
                    //situation that the particle is in the last column of the grid
                    else if (col == binNum-1){
                        for (int k = col-1; k <= col; k++){
                            for (int l = 0; l <bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                 
                            }
                        }
                    }
                }
            }

            //situation that the particle is in the last row of the grid
            if (row == binNum-1){
                for (int j = row-1; j <= row; j++){
                    //situation that the particle is not in the first or the last column of the grid
                    if ((col > 0) &&( col < binNum-1)){
                        for (int k = col-1; k <= col+1; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                  
                            }
                        }
                    }
                    //situation that the particle is in the first column of the grid
                    else if (col == 0){
                        for (int k = col; k <= col+1; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                 
                            }
                        }
                    }
                    //situation that the particle is in the last column of the grid
                    else{
                        for (int k = col-1; k <= col; k++){
                            for (int l = 0; l < bin[j*binNum + k].size(); l++){
                                 int fa = bin[j*binNum + k].at(l);
                                 apply_force(particles[i], particles[fa], &dmin, &davg, &navg);                                  
                            }
                        }
                    }
                }
            }
        }
        //
        //before moving the particles we need to release the inside vectors
        //
        for (int i = 0; i < binNum*binNum; i++)
            bin[i].resize(0);

        //
        //  move particles
        //

        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
