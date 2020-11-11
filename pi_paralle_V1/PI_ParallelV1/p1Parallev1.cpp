#include <iostream>
#include <unistd.h>
#include <chrono>
#include <omp.h>

static long num_steps=1000000;
double step;

#define NUM_THREADS 8


using namespace std;

int main()
{
    //double x;
    //double sum=0.0;
    //double pi;

    cout <<" steps ;" << ";"<< " Sum ;" << "PI ;" << "time" <<std::endl;
    const int num_testes = 50;
    for (int k= 0; k < num_testes; k++){


    int i, nthreads;
    double pi;
    double sum[NUM_THREADS]={0.0};
    step = 1.0/(double)num_steps;


     //Marcador inicial do tempo
    auto start_time = chrono::steady_clock::now();
    omp_set_num_threads(NUM_THREADS);

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


    pi=0.0;
    for (int i=0; i < nthreads; i++)
        pi += sum[i] *step ;

    auto end_time = chrono::steady_clock::now();

    cout << step << ";"<< sum << " ;" << pi << " ;"<< chrono::duration_cast<chrono::nanoseconds>(end_time-start_time).count() << endl;
    }
    return 0;
}
