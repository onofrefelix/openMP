#include <iostream>
#include <unistd.h>
#include <chrono>
#include <omp.h>

static long num_steps=100000000;
double step;

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




int main() {
    //double x;
    //double sum=0.0;
    //double pi;
     cout <<"PI ;" << "time Pal Crt ; " <<" PI ;" <<"Time  Seq" << std::endl;
     const int num_testes = 50;
     for (int k= 0; k < num_testes; k++){

      auto start_time_seq = chrono::steady_clock::now();
      double pi_s= calculaPISeq();
      auto end_time_seq = chrono::steady_clock::now();



    //Marcador inicial do tempo
        auto start_time = chrono::steady_clock::now();
        double pi = criticalParallel(num_steps,4);
        auto end_time = chrono::steady_clock::now();
        cout << pi << " ; "<< chrono::duration_cast<chrono::microseconds>(end_time-start_time).count() << " ; "<<
        pi << " ; "<< chrono::duration_cast<chrono::microseconds>(end_time_seq-start_time_seq).count() << endl;
   }
    return 0;
}
