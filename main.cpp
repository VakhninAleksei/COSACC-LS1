#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include "Header.h"
#include "lsgocec2013benchmarks.h"
#include "Constants.h"
using namespace std;

const int N = 1000, FEV_LS1_budget = 25000, FEV_global = 3e6, krantost = FEV_global/100, R = 3,
          generations_min = 5, generations_init = 20, pop_size_min = 25, pop_size_init = 100, pop_size_max = 200;
int ID, M, max_CC;
double a = -100.0, b = 100.0, best_solution = 1e300, THE_BEST_solution = 1e300, THE_GLOBAL_BEST_solution = 1e300;
const double LS1_step = 0.5,  piece = 0.1;

bool tournament_selection = 0;
const int H = 6;
int archive_size = pop_size_max*2, r1, r2, pbest, trigger_data_stat;
const int max_generations = 50;

string name_of_func;

int params_1[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
const int setups_number = 8;
int island_setup [setups_number][3] =
{
    {1,2,4}, {1,2,8}, {1,2,10}, {1,4,8}, {1,4,10}, {1,8,10}, {2,4,8}, {2,4,10}
};
int main(int argc, char** argv)
{
    int world_size, world_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Get_processor_name(processor_name, &name_len);


    int params_2[setups_number];

    for (int i=0; i!=setups_number; i++)
    {
        params_2[i]=i;
    }


    int params_3[1] = {0};
    int param_1_counter = sizeof(params_1)/sizeof(params_1[0]);
    int param_2_counter = sizeof(params_2)/sizeof(params_2[0]);

    int all_params = param_1_counter*param_2_counter;

    int *thread_number = new int [all_params];

    int j=0;

    for (int i=0; i!=all_params; i++)
    {
        thread_number[i] = j;
        j++;
        if (j > world_size-1)
        {
            j=0;
        }
    }

    int thread_index = 0;
    for (int p1=0; p1!=param_1_counter; p1++)
    {
        for (int p2=0; p2!=param_2_counter; p2++)
        {
            if (world_rank == thread_number[thread_index])
            {
                int island_pool[3];

                for (int i=0; i!=3; i++)
                {
                    island_pool[i] = island_setup[p2][i];
                }
                int pool_size=sizeof(island_pool)/sizeof(island_pool[0]);

                string settings = "{";
                for (int i=0; i!=pool_size; i++)
                {
                    if (i!=pool_size-1)
                    {
                        settings+=to_string(island_pool[i])+", ";
                    }
                    else
                    {
                        settings+=to_string(island_pool[i]);
                    }
                }

                settings +="}";
                cout<<settings<<endl;


                max_CC_parts (island_pool, max_CC, pool_size);

                int *FEV_island_generations  = new int [pool_size];
                int *FEV_island_generations_previous = new int [pool_size];

                srand(time(0));
                ID = params_1[p1];
                M = max_CC;

                cout<<"FEV: "<<FEV_global<<"| R: "<<R<< "| ID: "<<params_1[p1]<<"| M: "<<M<<"| pop_size: "<<pop_size_init<<" | settings: " <<settings<< " generations_init: " << generations_init<<" min_gen: "<<generations_min<<"  max_ger: "<<max_generations<<" tournament_selection: "<<tournament_selection<<endl<<endl;


                int **indices = new int*[pool_size];
                for (int count = 0; count < pool_size; count++)
                {
                    indices[count] = new int[N];
                }
                int S = N / M;
                int *range = new int[M + 1];
                range[0] = 0;
                range[M] = N;
                for (int i = 1; i < M; i++)
                {
                    range[i] = range[i - 1] + S;
                }

                double *Ovector = new double[N];
                int *Pvector = new int[N];

                double **r25 = new double*[25];
                for (int count = 0; count < 25; count++)
                {
                    r25[count] = new double[25];
                }

                double **r50 = new double*[50];
                for (int count = 0; count < 50; count++)
                {
                    r50[count] = new double[50];
                }

                double **r100 = new double*[100];
                for (int count = 0; count < 100; count++)
                {
                    r100[count] = new double[100];
                }

                int s_size = 0;

                if (ID == 4 || ID == 5 || ID == 6 || ID == 7)
                {
                    s_size = 7;
                }
                if (ID == 8 || ID == 9 || ID == 10 || ID == 11 || ID == 13 || ID == 14)
                {
                    s_size = 20;
                }

                int *s = new int[s_size];
                double *w = new double[s_size];

                if (ID == 4 || ID == 5 || ID == 6 || ID == 7 || ID == 8 || ID == 9 || ID == 10 || ID == 11 || ID == 13 || ID == 14)
                {
                    Pvector = readPermVector(N, ID);
                    r25 = readR(25, ID);
                    r50 = readR(50, ID);
                    r100 = readR(100, ID);
                    s = readS(s_size, ID);
                    w = readW(s_size, ID);
                }
                double **OvectorVec = new double*[s_size];
                for (int count = 0; count < s_size; count++)
                    OvectorVec[count] = new double[s[count]];

                if (ID == 14)
                {
                    OvectorVec = readOvectorVec(N, s_size, s, ID);
                }
                if (ID != 14)
                {
                    Ovector = readOvector(N, ID);
                }

                select_borders(a, b, name_of_func, ID);

                string piece_string = to_string(piece);
                piece_string.erase ( piece_string.find_last_not_of('0') + 1, std::string::npos );
                int *k = new int [M];
                int **A = new int*[pool_size];
                for (int count = 0; count < pool_size; count++)
                {
                    A[count] = new int[M];
                }

                string name ="Experiment_results/COSACC_LS1_"+ to_string(ID) +"_"+settings
                +"_("+to_string(generations_min)+"_"+to_string(generations_init)+")_"+to_string(pop_size_init)+
                "_"+to_string(pop_size_min)+"_"+to_string(pop_size_max)+".txt";
                ofstream fout(name);
                fout << "COSACC-LS1 LSGO CEC'2013"<< endl<< "ID: " << ID << endl<<"FEV: "<<FEV_global<<endl<<"Name_of benchmark problem: "<<name_of_func<<endl;
                fout << "R: "<< R << endl <<"Population size: "<< pop_size_init<<endl<< "N: " << N <<endl << "a: " << a << endl << "b: " << b << endl;
                fout << "M: " << M <<endl<<"H: "<<H<<endl<<"Archive_size: "<<archive_size<<endl<<"Piece: "<<piece<<" settings: "<<settings<<" generations_init: "<<generations_init<<" tournament_selection: "<<tournament_selection<<" LS1: "<<FEV_LS1_budget<<endl;


                double **population = new double *[pop_size_max];
                for(int j=0; j<pop_size_max; j++)
                {
                    population[j] = new double [N];
                }

                double **population_new = new double *[pop_size_max];
                for(int j=0; j<pop_size_max; j++)
                {
                    population_new[j] = new double [N];
                }

                double ***fitness_cc= new double **[pool_size];
                for(int i=0; i<pool_size; i++)
                {
                    fitness_cc[i] = new double *[M];
                    for(int j=0; j<M; j++)
                    {
                        fitness_cc[i][j] = new double [pop_size_max];
                    }
                }

                double ***FEV_island_stat= new double **[R];
                for(int i=0; i<R; i++)
                {
                    FEV_island_stat[i] = new double *[pool_size];
                    for(int j=0; j<pool_size; j++)
                    {
                        FEV_island_stat[i][j] = new double [200];
                    }
                }

                int **pop_size_island_stat= new int *[R];
                for(int i=0; i<R; i++)
                {
                    pop_size_island_stat[i] = new int [200];
                }

                double ***fitness_cc_new= new double **[pool_size];
                for(int i=0; i<pool_size; i++)
                {
                    fitness_cc_new[i] = new double *[M];
                    for(int j=0; j<M; j++)
                    {
                        fitness_cc_new[i][j] = new double [pop_size_max];
                    }
                }

                double **archive = new double *[archive_size];
                for(int j=0; j<archive_size; j++)
                {
                    archive[j] = new double [N];
                }

                double ***HISTORY_F= new double **[pool_size];
                for(int i=0; i<pool_size; i++)
                {
                    HISTORY_F[i] = new double *[M];
                    for(int j=0; j<M; j++)
                    {
                        HISTORY_F[i][j] = new double [H];
                    }
                }
                double ***HISTORY_CR= new double **[pool_size];
                for(int i=0; i<pool_size; i++)
                {
                    HISTORY_CR[i] = new double *[M];
                    for(int j=0; j<M; j++)
                    {
                        HISTORY_CR[i][j] = new double [H];
                    }
                }

                int **cc_best_individual_index = new int *[pool_size];
                for (int count = 0; count < pool_size; count++)
                {
                    cc_best_individual_index[count] = new int[M];
                }


                double *S_CR = new double [pop_size_max];
                double *S_F = new double [pop_size_max];
                double *delta_f = new double [pop_size_max];
                double *W = new double [pop_size_max];
                double *solution = new double [N];
                double *solution_copy = new double [N];
                int *r = new int [pop_size_max];
                double *F = new double [pop_size_max];
                double *CR = new double [pop_size_max];

                double **u = new double *[pop_size_max];
                for (int count = 0; count < pop_size_max; count++)
                {
                    u[count] = new double [N];
                }

                double **data_stat = new double *[R];
                for (int count = 0; count < R; count++)
                {
                    data_stat[count] = new double [200];
                }

                double **best_solution_ever = new double *[R];
                for (int count =0; count <R; count++)
                {
                    best_solution_ever[count] = new double [N];
                }
                int best_solution_ever_M;
                int *best_solution_ever_CC_individuals = new int[max_CC];

                double *best_fitness_ever = new double [R];
                for (int i=0; i!=R; i++)
                {
                    best_fitness_ever[i] = 1e300;
                }
                int **data_restart_stat = new int *[R];
                for (int count = 0; count < R; count++)
                {
                    data_restart_stat[count] = new int [2000];
                }
                int *restart_trigger = new int [R];
                for(int i=0; i!=R; i++)
                {
                    restart_trigger[i]=0;
                }

                int *sequense_of_islands = new int [pool_size];
                int *extra_pool = new int [pool_size];
                for (int i=0; i!=pool_size; i++)
                {
                    sequense_of_islands[i]=i;
                }

                for (int i=0; i!=pool_size; i++)
                {
                    extra_pool[i] = island_pool[i];
                }


                for (int z=0; z!=R; z++)
                {
                    double *SR = new double[N];
                    for (int j=0; j!=N; j++)
                    {
                        SR[j] = 0.4*(b-a);
                    }
                    int *LS1_indecies = new int [N];
                    double *LS1_improve = new double [N];
                    indecesSuccession(LS1_indecies, N);
                    randperm(LS1_indecies, N);

                    THE_BEST_solution = 1e300;
                    THE_GLOBAL_BEST_solution = 1e300;

                    int restarts = 0;
                    bool trigger_restart = 0;
                    int num_worse = 0;
                    int pop_size = pop_size_init;
                    indecesSuccession(indices, pool_size, N);
                    randperm_in_all_islands(indices, N, pool_size);

                    filling_zeros_A(A, M, pool_size);
                    filling_zeros(k, M);
                    filling_zeros(S_CR, pop_size_max);
                    filling_zeros(S_F, pop_size_max);
                    filling_zeros(delta_f, pop_size_max);

                    trigger_data_stat = 0;

                    int FEV = FEV_global;
                    int success = 0;

                    for (int i=0; i!=pool_size; i++)
                    {
                        pop_size = pop_size_init;
                    }

                    initializePopulation(population, population_new, pop_size, N, a, b);
                    double DI_init = DI_calculation(population, pop_size, N);

                    bubble_sort_indecies(extra_pool, sequense_of_islands, pool_size);
                    for (int island=0; island!=pool_size; island++)
                    {
                        M = island_pool[island];
                        initializeHistory(HISTORY_F, HISTORY_CR, H, M, island);
                    }

                    for (int i=0; i!=pool_size; i++)
                    {
                        FEV_island_generations[i] = FEV_island_generations_previous [i] = generations_init;
                    }

                    double *previous_best_fitness_island = new double [pool_size];
                    double *current_best_fitness_island = new double [pool_size];
                    double *island_performance = new double [pool_size];
                    double previous_Improvement;
                    double current_Improvement;
                    double best_for_LS1 = 1e300;

                    while (FEV>0)
                    {
                        bool improvement_trigger = 0;
                        previous_Improvement = 1e300;
                        current_Improvement = 1e300;
                        randperm(sequense_of_islands, pool_size);
                        for (int island_original=0; island_original!=pool_size; island_original++)
                        {
                            int island = sequense_of_islands[island_original];
                            M = island_pool[island];

                            THE_BEST_solution = 1e300;
                            rnd_indecies(cc_best_individual_index, pop_size, M, island);

                            recalculating_range_interval(range, M, N, S);

                            for (int p=0; p!=M; p++)
                            {
                                for (int i=0; i!=pop_size; i++)
                                {
                                    building_solution(range, solution, indices, cc_best_individual_index, population, i, p, M, island);
                                    for (int k=0; k!=N; k++)
                                    {
                                        solution_copy[k] = solution[k];
                                    }
                                    fitness_cc[island][p][i] = fitness_cc_new[island][p][i] = benchmark_func(solution, Ovector, OvectorVec, Pvector, r25, r50, r100, s, w, N, ID);
                                    if (fitness_cc[island][p][i]<best_for_LS1)
                                    {
                                        best_for_LS1 = fitness_cc[island][p][i];
                                        for (int k=0; k!=N; k++)
                                        {
                                            best_solution_ever[z][k] = solution_copy[k];
                                        }
                                    }
                                    FEV--;
                                    if (THE_GLOBAL_BEST_solution > fitness_cc[island][p][i])
                                    {
                                        THE_GLOBAL_BEST_solution = fitness_cc[island][p][i];
                                    }
                                    if (THE_BEST_solution > fitness_cc[island][p][i])
                                    {
                                        THE_BEST_solution = fitness_cc[island][p][i];
                                    }

                                    if (FEV % krantost == 0 )
                                    {

                                        data_stat[z][trigger_data_stat] = THE_GLOBAL_BEST_solution;
                                        best_fitness_ever[z] = THE_GLOBAL_BEST_solution;

                                        double sum = 0.0;
                                        for (int k=0; k!=pool_size; k++)
                                        {
                                            sum+=FEV_island_generations_previous[k];
                                        }
                                        for (int k=0; k!=pool_size; k++)
                                        {
                                            FEV_island_stat[z][k][trigger_data_stat] = (double)FEV_island_generations_previous[k]/sum;
                                            pop_size_island_stat[z][trigger_data_stat] = pop_size;
                                        }
                                        trigger_data_stat++;
                                    }
                                }
                                find_best_part_index(cc_best_individual_index,fitness_cc, p, pop_size, island);
                            }

                            if (improvement_trigger == 0)
                            {
                                improvement_trigger = 1;
                                previous_Improvement = median_island(1, pop_size, island, M, fitness_cc);
                                current_Improvement = current_Improvement;
                            }

                            previous_best_fitness_island[island] = median_island(1, pop_size, island, M, fitness_cc);

                            while (FEV_island_generations[island]>0)
                            {
                                for (int p=0; p!=M; p++)
                                {
                                    for (int i=0; i!=pop_size; i++)
                                    {
                                        GENERATING_TRIAL_VECTOR (i, fitness_cc, pbest, pop_size, piece, p, r, H, F, CR, HISTORY_CR, HISTORY_F, r1, r2, A, range, u, population, indices, archive, island, tournament_selection);
                                        check_out_borders(u, population, i, N, a, b, range, p, indices, island);
                                    }

                                    for (int i=0; i!=pop_size; i++)
                                    {
                                        building_solution_U(range, solution, u, indices, cc_best_individual_index, population, i, p, M, island);
                                        for (int k=0; k!=N; k++)
                                        {
                                            solution_copy[k] = solution[k];
                                        }
                                        double test_function = benchmark_func(solution, Ovector, OvectorVec, Pvector, r25, r50, r100, s, w, N, ID);
                                        FEV--;
                                        if (test_function < best_for_LS1)
                                        {
                                            best_for_LS1 = test_function;
                                            for (int k=0; k!=N; k++)
                                            {
                                                best_solution_ever[z][k] = solution_copy[k];
                                            }
                                        }
                                        if (current_Improvement > test_function)
                                        {
                                            current_Improvement = test_function;
                                        }

                                        if (THE_GLOBAL_BEST_solution > test_function)
                                        {
                                            THE_GLOBAL_BEST_solution = test_function;
                                            best_fitness_ever[z] = THE_GLOBAL_BEST_solution;
                                        }
                                        if (THE_BEST_solution > test_function)
                                        {
                                            THE_BEST_solution = test_function;
                                        }

                                        if (current_best_fitness_island[island] > test_function)
                                        {
                                            current_best_fitness_island[island] = test_function;
                                        }

                                        if (fitness_cc[island][p][i] > test_function)
                                        {
                                            for (int j=range[p]; j!=range[p+1]; j++)
                                            {
                                                population_new[i][indices[island][j]] = u[i][indices[island][j]];
                                            }
                                            UPDATING_INFO(archive, population, i, archive_size, A, range, p, indices, fitness_cc, fitness_cc_new, test_function, delta_f, S_F, F, S_CR, CR, success, island);
                                        }

                                        if (FEV % krantost == 0 )
                                        {
                                            information_output (name_of_func, THE_GLOBAL_BEST_solution,THE_BEST_solution,FEV, z, R, M, pop_size, archive_size, H, HISTORY_CR, HISTORY_F, island, restarts);
                                            data_stat[z][trigger_data_stat] = THE_GLOBAL_BEST_solution;
                                            best_fitness_ever[z] = THE_GLOBAL_BEST_solution;
                                            double sum = 0.0;
                                            for (int k=0; k!=pool_size; k++)
                                            {
                                                sum+=FEV_island_generations_previous[k];
                                            }
                                            for (int k=0; k!=pool_size; k++)
                                            {
                                                FEV_island_stat[z][k][trigger_data_stat] = (double)FEV_island_generations_previous[k]/sum;
                                                pop_size_island_stat[z][trigger_data_stat] = pop_size;
                                            }

                                            trigger_data_stat++;
                                        }
                                    }

                                    Algorithm_1(delta_f, W, S_CR, S_F, HISTORY_CR, HISTORY_F, k, success, H, p, island);

                                    for (int i=0; i!=pop_size; i++)
                                    {
                                        for (int j=range[p]; j!=range[p+1]; j++)
                                        {
                                            population[i][indices[island][j]] = population_new[i][indices[island][j]];
                                        }
                                        fitness_cc[island][p][i] = fitness_cc_new[island][p][i];
                                    }
                                    find_best_part_index(cc_best_individual_index, fitness_cc, p, pop_size, island);
                                }

                                FEV_island_generations[island]--;
                                current_best_fitness_island[island] = median_island(0.05, pop_size, island, M, fitness_cc);

                            }

                            if (previous_best_fitness_island[island] == 0)
                            {
                                island_performance[island] = 1.0;
                            }
                            if (previous_best_fitness_island[island] > 0)
                            {
                                island_performance[island] = (previous_best_fitness_island[island]-current_best_fitness_island[island])/previous_best_fitness_island[island];
                                island_performance[island] /=FEV_island_generations_previous[island];
                            }


                            DI_Algorithm(population, population_new, pop_size, N, DI_init, FEV_global, FEV, pop_size_min, pop_size_max,  a, b, island_pool, range, S, fitness_cc, fitness_cc_new, island);


                        }

                        current_Improvement = median_all_islands (1, pop_size, fitness_cc, pool_size, island_pool);

                        Improvement(previous_Improvement, current_Improvement, 0.05, num_worse, trigger_restart);
                        FEV_island_calculating(FEV_island_generations, FEV_island_generations_previous, generations_min, max_generations, island_performance, pool_size);

                        for (int i=0; i!=pool_size; i++)
                        {
                            cout<<i<<": "<<island_performance[i]<<" "<<FEV_island_generations[i]<<endl;
                        }




                        int FEV_LS1 = FEV_LS1_budget;
                        bool best_improved = 0;
                        double improvement_LS1 = best_for_LS1;
                        double LS1_fitness_before = best_for_LS1;
                        filling_zeros(LS1_improve, N);

                        while (FEV_LS1 > 0)
                        {
                            filling_zeros(LS1_improve, N);
                            for (int i=0; i!=N; i++)
                            {
                                double current_best_fitness = best_for_LS1;
                                double *newsol = new double [N];
                                for (int j=0; j!=N; j++)
                                {
                                    newsol[j] = best_solution_ever[z][j];
                                }
                                newsol[LS1_indecies[i]] -= SR[LS1_indecies[i]];
                                borbers_check(newsol, N, a, b);

                                double test_function = benchmark_func(newsol, Ovector, OvectorVec, Pvector, r25, r50, r100, s, w, N, ID);
                                FEV--;
                                FEV_LS1--;

                                if (test_function < THE_GLOBAL_BEST_solution)
                                {
                                    THE_GLOBAL_BEST_solution = test_function;
                                    best_fitness_ever[z] = THE_GLOBAL_BEST_solution;

                                }
                                if (test_function < best_for_LS1)
                                {
                                    best_for_LS1 = test_function;
                                    for (int q=0; q!=N; q++)
                                    {
                                        best_solution_ever[z][q] = newsol[q];
                                    }
                                    best_improved = 1;

                                }

                                if (FEV % krantost == 0 )
                                {

                                    data_stat[z][trigger_data_stat] = THE_GLOBAL_BEST_solution;
                                    best_fitness_ever[z] = THE_GLOBAL_BEST_solution;

                                    double sum = 0.0;
                                    for (int k=0; k!=pool_size; k++)
                                    {
                                        sum+=FEV_island_generations_previous[k];
                                    }
                                    for (int k=0; k!=pool_size; k++)
                                    {
                                        FEV_island_stat[z][k][trigger_data_stat] = (double)FEV_island_generations_previous[k]/sum;
                                        pop_size_island_stat[z][trigger_data_stat] = pop_size;
                                    }
                                    trigger_data_stat++;
                                }

                                if (test_function > best_for_LS1)
                                {
                                    for (int j=0; j!=N; j++)
                                    {
                                        newsol[j] = best_solution_ever[z][j];
                                    }

                                    newsol[LS1_indecies[i]] +=LS1_step*SR[LS1_indecies[i]];

                                    borbers_check(newsol, N, a, b);
                                    test_function = benchmark_func(newsol, Ovector, OvectorVec, Pvector, r25, r50, r100, s, w, N, ID);
                                    FEV--;
                                    FEV_LS1--;

                                    if (test_function < THE_GLOBAL_BEST_solution)
                                    {
                                        THE_GLOBAL_BEST_solution = test_function;
                                        best_fitness_ever[z] = THE_GLOBAL_BEST_solution;
                                    }

                                    if (test_function < best_for_LS1)
                                    {
                                        best_for_LS1= test_function;
                                        for (int q=0; q!=N; q++)
                                        {
                                            best_solution_ever[z][q] = newsol[q];
                                        }
                                        best_improved = 1;

                                    }

                                    if (FEV % krantost == 0 )
                                    {
                                        data_stat[z][trigger_data_stat] = THE_GLOBAL_BEST_solution;
                                        best_fitness_ever[z] = THE_GLOBAL_BEST_solution;

                                        double sum = 0.0;
                                        for (int k=0; k!=pool_size; k++)
                                        {
                                            sum+=FEV_island_generations_previous[k];
                                        }
                                        for (int k=0; k!=pool_size; k++)
                                        {
                                            FEV_island_stat[z][k][trigger_data_stat] = (double)FEV_island_generations_previous[k]/sum;
                                            pop_size_island_stat[z][trigger_data_stat] = pop_size;
                                        }

                                        trigger_data_stat++;
                                    }
                                }
                                delete [] newsol;
                                LS1_improve[LS1_indecies[i]] = std::max(improvement_LS1 - best_for_LS1, 0.0);
                                if (improvement_LS1>best_for_LS1)
                                {
                                    improvement_LS1 = best_for_LS1;
                                }
                                LS_SR_update(i, LS1_indecies, LS1_improve, SR, a, b);
                            }
                            bubble_sort_indecies(LS1_improve, LS1_indecies, N);
                        }
                        double LS1_fitness_after = best_for_LS1;

                        if (LS1_fitness_after == 0.0)
                        {

                        }
                        else
                        {
                            Improvement(LS1_fitness_before, LS1_fitness_after, 0.05, num_worse, trigger_restart);
                        }

                        if (best_improved == 1)
                        {
                            int rnd = RANDOM()*(pop_size-1);

                            for (int j=0; j!=N; j++)
                            {
                                population[rnd][j] = population_new[rnd][j] = best_solution_ever[z][j];
                            }

                        }
                        else
                        {
                            randperm(LS1_indecies, N);
                        }


                        best_for_LS1 = 1e300;

                        double aver_SR = 0;
                        for (int i=0; i!=N; i++)
                        {
                            aver_SR +=SR[i];
                        }


                        if (LS1_fitness_before > 0)
                        {
                            double imp = (LS1_fitness_before-LS1_fitness_after)/LS1_fitness_before;
                            if (imp <0.05)
                            {
                                for (int j=0; j!=N; j++)
                                {
                                    SR[j] = 0.4*(b-a);
                                }
                            }

                        }




                    }

                }

                fout<<endl;

                for (int i=0; i!=R; i++)
                {
                    for (int j=0; j!=100; j++)
                    {
                        if (j<99)
                        {
                            fout<<data_stat[i][j]<<", ";
                        }

                        if (j == 99)
                        {
                            fout<<data_stat[i][j];
                        }
                    }
                    fout<<endl;
                }

                fout << endl;
                fout << "Average:" << endl;
                for (int i = 0; i != 100; i++)
                {
                    double sum = 0.0;

                    for (int j = 0; j != R; j++)
                    {
                        sum += data_stat[j][i];
                    }

                    if (i<99)
                    {
                        fout << sum / R << ", ";
                    }

                    if (i==99)
                    {
                        fout << sum / R;
                    }
                }

                fout << endl;
                fout << endl;
                fout << "Best decisions" << endl;

                for (int i = 0; i < R; i++)
                {
                    data_stat[i][99] = best_fitness_ever[i];
                    if (i!=(R-1))
                    {
                        fout << data_stat[i][99] << ", ";
                    }
                    if (i == R-1)
                    {
                        fout << data_stat[i][99];
                    }
                }
                fout << endl;
                fout << endl;
                int sr1 = 4;
                fout << endl;
                fout << "1.2e+5: " << endl;
                fout << "BEST: " << min_stat(data_stat, sr1 - 1, R) << endl;
                fout << "MEDIAN: " << median_stat(data_stat, sr1 - 1, R) << endl;
                fout << "WORST: " << max_stat(data_stat, sr1 - 1, R) << endl;
                fout << "MEAN: " << mean_stat(data_stat, sr1 - 1, R) << endl;
                fout << "STDEV: " << stddev_stat(data_stat, sr1 - 1, R, mean_stat(data_stat, sr1 - 1, R)) << endl;
                fout << endl;

                int sr2 = 20;
                fout << endl;
                fout << "6.0e+5: " << endl;
                fout << "BEST: " << min_stat(data_stat, sr2 - 1, R) << endl;
                fout << "MEDIAN: " << median_stat(data_stat, sr2 - 1, R) << endl;
                fout << "WORST: " << max_stat(data_stat, sr2 - 1, R) << endl;
                fout << "MEAN: " << mean_stat(data_stat, sr2 - 1, R) << endl;
                fout << "STDEV: " << stddev_stat(data_stat, sr2 - 1, R, mean_stat(data_stat, sr2 - 1, R)) << endl;
                fout << endl;


                int sr3 = 100;
                fout << endl;
                fout << "3e+6: " << endl;
                fout << "BEST: " << min_stat(data_stat, sr3 - 1, R) << endl;
                fout << "MEDIAN: " << median_stat(data_stat, sr3 - 1, R) << endl;
                fout << "WORST: " << max_stat(data_stat, sr3 - 1, R) << endl;
                fout << "MEAN: " << mean_stat(data_stat, sr3 - 1, R) << endl;
                fout << "STDEV: " << stddev_stat(data_stat, sr3 - 1, R, mean_stat(data_stat, sr3 - 1, R)) << endl;
                fout << endl;
                fout << endl;

                fout << endl;
                fout << min_stat(data_stat, sr3 - 1, R) << endl;
                fout << median_stat(data_stat, sr3 - 1, R) << endl;
                fout << max_stat(data_stat, sr3 - 1, R) << endl;
                fout << mean_stat(data_stat, sr3 - 1, R) << endl;
                fout << stddev_stat(data_stat, sr3 - 1, R, mean_stat(data_stat, sr3 - 1, R)) << endl;

                for (int i=0; i!=R; i++)
                {
                    fout<<i+1<<": "<<endl;
                    for (int k=0; k!=pool_size; k++)
                    {
                        for (int j=0; j!=100; j++)
                        {
                            if (j<99)
                            {
                                fout<<FEV_island_stat[i][k][j]<<", ";
                            }

                            if (j == 99)
                            {
                                fout<<FEV_island_stat[i][k][j];
                            }
                        }
                        fout<<endl;
                    }

                    fout<<endl;
                    fout<<endl;
                }

                fout<<"AVERAGE: "<<endl;

                for (int j=0; j!=pool_size; j++)
                {
                    for (int i=0; i!=100; i++)
                    {
                        double SUM = 0.0;
                        if (i<99)
                        {


                            for (int k=0; k!=R; k++)
                            {
                                SUM+=FEV_island_stat[k][j][i];
                            }
                            fout<<SUM/R<<", ";
                        }
                        if (i==99)
                        {


                            for (int k=0; k!=R; k++)
                            {
                                SUM+=FEV_island_stat[k][j][i];
                            }
                            fout<<SUM/R;
                        }
                    }
                    fout<<endl;

                }

                for (int i=0; i!=R; i++)
                {
                    fout<<i+1<<": "<<endl;

                    for (int j=0; j!=100; j++)
                    {
                        if (j<99)
                        {
                            fout<<pop_size_island_stat[i][j]<<", ";
                        }

                        if (j == 99)
                        {
                            fout<<pop_size_island_stat[i][j];
                        }
                    }
                    fout<<endl;
                    fout<<endl;
                    fout<<endl;
                }

                fout<<"AVERAGE_pop_size: "<<endl;


                for (int i=0; i!=100; i++)
                {

                    double SUM = 0.0;
                    if (i<99)
                    {
                        for (int k=0; k!=R; k++)
                        {
                            SUM+=pop_size_island_stat[k][i];
                        }
                        fout<<SUM/R<<", ";
                    }

                    if (i==99)
                    {
                        for (int k=0; k!=R; k++)
                        {
                            SUM+=pop_size_island_stat[k][i];
                        }
                        fout<<SUM/R;
                    }

                }
                fout<<endl;

                fout.close();



            }
            thread_index++;

        }
    }
    MPI_Finalize();
}
