#include <iostream>

#include <unistd.h>

#include <chrono>

#include <omp.h>



static long num_steps=1000000;



double step;



using namespace std;



int main()

{

    //omp_set_num_threads(1);



    double x;

    double sum=0.0, locsum=0.0;

    double pi;

    int numthreads;





     //Marcador inicial do tempo

    auto start_time = chrono::steady_clock::now();

    step = 1.0/(double)num_steps;



    omp_set_num_threads(4);



    numthreads = omp_get_num_threads();

    cout << numthreads << "\n";



    for(int k = 0; k < 10; k++){

    sum = 0; locsum = 0;



        #pragma omp parallel private(x, locsum)

        {

            numthreads = omp_get_num_threads();



            #pragma omp for

            for (int i=0; i < num_steps; i++)

            {

                x = (i+0.5)*step;

                locsum += 4.0/(1.0+x*x);

            }



            #pragma omp critical

            {

                sum += locsum;

            }

        }

        pi = step * sum;



        auto end_time = chrono::steady_clock::now();



        cout << "threads " << numthreads <<" steps = " << step << " Sum = "<< sum << " PI = " << pi << "time = "<< chrono::duration_cast<chrono::nanoseconds>(end_time-start_time).count() << " ns" << endl;



    }

    return 0;



}
