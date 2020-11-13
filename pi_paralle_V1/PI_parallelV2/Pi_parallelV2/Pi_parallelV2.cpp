#include <iostream>
#include <unistd.h>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <climits>

static long num_steps=10000000;
double step;

const double PI=3.141592653589793238;

#define NUM_THREADS 4


using namespace std;

double calculaPISeq()
{

    double x;
    double sum=0.0;
    double pi;
     //Marcador inicial do tempo
    step = 1.0/(double)num_steps;

    for (int i=0; i < num_steps; i++)
    {
        x = (i+0.5)*step;
        sum += 4.0/(1.0+x*x);
    }

    pi = step * sum;
    return pi;

}



double reduction(long int num_steps, int numThreads)
{
     int i; double x, pi, sum = 0.0;
     step = 1.0/(double) num_steps;
     omp_set_num_threads(numThreads);
     #pragma omp parallel for private(x) reduction(+:sum)
     for (i=0;i< num_steps; i++){
        x = (i+0.5)*step;
        sum = sum + 4.0/(1.0+x*x);
        }
     pi = step * sum;

     return pi;

}


double falseSharing(long int num_steps, int numThreads)
{

    int i, nthreads;
    double pi=0.0;
    double sum[numThreads]={0.0};
    step = 1.0/(double)num_steps;

     omp_set_num_threads(numThreads);

    #pragma omp parallel
    {
        int i, id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if(id == 0) nthreads = nthrds;
        for (int i=id; i < num_steps; i+=nthrds)
        {
            x = (i+0.5)*step;
            sum[id] += 4.0/(1.0+x*x);
        }
    }


    for (int i=0; i < nthreads; i++)
        pi += sum[i] *step ;

    return pi;
}





double criticalParallel(long int num_steps, int numThreads)

{
   omp_set_num_threads(numThreads);
   int i, nthreads;
   double pi=0.0;

   step = 1.0/(double)num_steps;


    #pragma omp parallel
    {
        int i, id,nthrds;
        double x, sum=0.0;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if(id == 0) nthreads = nthrds;

        for (int i=id; i < num_steps; i+=nthrds)
        {
            x = (i+0.5)*step;
            sum += 4.0/(1.0+x*x);
        }
        #pragma omp critical
            pi += sum*step;
    }


    return pi;
}


double atomic_Parallel(long int num_steps, int numThreads)

{
   omp_set_num_threads(numThreads);
   int i, nthreads;
   double pi=0.0;

   step = 1.0/(double)num_steps;


    #pragma omp parallel
    {
        int i, id,nthrds;
        double x, sum=0.0;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if(id == 0) nthreads = nthrds;

        for (int i=id; i < num_steps; i+=nthrds)
        {
            x = (i+0.5)*step;
            sum += 4.0/(1.0+x*x);
        }
        #pragma omp atomic
            pi += sum*step;
    }

    return pi;
}





int main() {
    //double x;
    //double sum=0.0;
    //double pi;
     cout << " num_threads ; " <<" PI_s ;" << " Relative Error ;"  << "Time  serial ;" << "speed up ;"
                               <<"PI_re ;" << " Relative Error ;"  << "time Pal re ; " << "speed up ;"
                               <<"PI_fs ;" << " Relative Error ;"  << "time Pal fs; " << "speed up ;"
                               <<"PI_a  ;" << " Relative Error ;"  << "time Pal atom ; " << "speed up ;"
                               <<"PI_c ;"  << " Relative Error ;"  << "Time  crit " << "speed up" << std::endl;
     const int num_testes = 20;
     double ttime_re = 0.0;
     double ttime_s  = 0.0;
     double ttime_c  = 0.0;
     double ttime_fs = 0.0;
     double ttime_a  = 0.0;
     double pi_s = 0.0;
     double pi_re = 0.0;
     double pi_fs = 0.0;
     double pi_c = 0.0;
     double pi_a = 0.0;


     for (int numtrhead = 0; numtrhead < 16; numtrhead++ ){

         for (int k= 0; k < num_testes; k++){


            // calculo do PI com reduction
            auto start_time_re = chrono::steady_clock::now();
            pi_re += reduction(num_steps,numtrhead+1 );
            auto end_time_re = chrono::steady_clock::now();

            ttime_re += chrono::duration_cast<chrono::microseconds>(end_time_re-start_time_re).count();



              /// calculo do PI serial
             auto start_time_s = chrono::steady_clock::now();
             pi_s += calculaPISeq();
             auto end_time_s = chrono::steady_clock::now();

             ttime_s += chrono::duration_cast<chrono::microseconds>(end_time_s-start_time_s).count();


            /// calculo do PI com  critical
             auto start_time_c = chrono::steady_clock::now();
             pi_c += criticalParallel(num_steps,numtrhead+1 );
             auto end_time_c = chrono::steady_clock::now();

             ttime_c += chrono::duration_cast<chrono::microseconds>(end_time_c-start_time_c).count();

           /// calculo do PI com  Atomic
             auto start_time_a = chrono::steady_clock::now();
             pi_a += atomic_Parallel(num_steps,numtrhead+1);
             auto end_time_a = chrono::steady_clock::now();

             ttime_a += chrono::duration_cast<chrono::microseconds>(end_time_a-start_time_a).count();

                /// calculo do PI false sharing
             auto start_time_fs = chrono::steady_clock::now();
             pi_fs += falseSharing(num_steps,numtrhead+1);
             auto end_time_fs = chrono::steady_clock::now();

             ttime_fs += chrono::duration_cast<chrono::microseconds>(end_time_fs-start_time_fs).count();


               /*

                cout << numtrhead+1 <<" ;" << pi_s <<" ; " << fabs(pi_s - PI)/fabs(PI) << " ; "<< chrono::duration_cast<chrono::microseconds>(end_time_s-start_time_s).count() << " ; "
                     << pi_re <<" ; " << fabs(pi_re - PI)/fabs(PI) << " ; "<< chrono::duration_cast<chrono::microseconds>(end_time_re-start_time_re).count() << " ; "
                     << pi_fs <<" ; " << fabs(pi_fs - PI)/fabs(PI) << " ; "<< chrono::duration_cast<chrono::microseconds>(end_time_fs-start_time_fs).count() << " ; "
                     << pi_a <<" ; " << fabs(pi_a - PI)/fabs(PI) << " ; "<< chrono::duration_cast<chrono::microseconds>(end_time_a-start_time_a).count() << " ; "
                     << pi_c <<" ; " << fabs(pi_c - PI)/fabs(PI)<< " ; "<< chrono::duration_cast<chrono::microseconds>(end_time_c-start_time_c).count() << endl;  */
        }

        cout << numtrhead+1 <<" ;" << pi_s/num_testes <<" ; " << fabs(pi_s/num_testes - PI)/fabs(PI) << " ; "<< ttime_s/num_testes << " ; " <<  (ttime_s)/(ttime_s) << " ; "
                     << pi_re/num_testes <<" ; " << fabs(pi_re/num_testes - PI)/fabs(PI) << " ; "<< ttime_re/num_testes << " ; "   <<  (ttime_s)/(ttime_re) << " ; "
                     << pi_fs/num_testes <<" ; " << fabs(pi_fs/num_testes - PI)/fabs(PI) << " ; "<< ttime_fs/num_testes << " ; "   <<  (ttime_s)/(ttime_fs) << " ; "
                     << pi_a/num_testes <<" ; " << fabs(pi_a/num_testes - PI)/fabs(PI) << " ; "<< ttime_a/num_testes    << " ; "   <<  (ttime_s)/(ttime_a) << " ; "
                     << pi_c/num_testes <<" ; " << fabs(pi_c/num_testes - PI)/fabs(PI)<< " ; "<< ttime_c/num_testes     << " ; "   <<  (ttime_s)/(ttime_c) << " ; " <<  endl;

        pi_s = 0.0;
        pi_re = 0.0;
        pi_fs = 0.0;
        pi_c = 0.0;
        pi_a = 0.0;


     }
    return 0;
}
