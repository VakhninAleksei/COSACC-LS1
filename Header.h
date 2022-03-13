//#ifndef HEADER_H_INCLUDED
//#define HEADER_H_INCLUDED

#include <iostream>
#include <random>
#include <cstdlib>
#include <math.h>
#include "Constants.h"


using namespace std;

void quickSort(double *arr, int left, int right)
{

    int i = left, j = right;
    double tmp;
    double pivot = arr[(left + right) / (int)2];
    while (i <= j)
    {

        while (arr[i] < pivot)
            i++;

        while (arr[j] > pivot)
            j--;

        if (i <= j)
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    };

    if (left < j)
        quickSort(arr, left, j);

    if (i < right)
        quickSort(arr, i, right);

}

void bubble_sort(double *a, int lenght)
{
    for (int j=0; j<lenght -1; j++)
    {
        for (int i=0; i<lenght-j-1; i++)
        {
            if (a[i] > a[i+1])
            {
                double b = a[i];
                a[i] = a[i+1];
                a[i+1]=b;
            }
        }
    }
}
void bubble_sort_inverse(double *a, int lenght)
{
    for (int j=0; j<lenght -1; j++)
    {
        for (int i=0; i<lenght-j-1; i++)
        {
            if (a[i] < a[i+1])
            {
                double b = a[i];
                a[i] = a[i+1];
                a[i+1]=b;
            }
        }
    }
}

void bubble_sort_indecies(double *a, int *index, int lenght)
{
    for (int j=0; j<lenght-1; j++)
    {
        for (int i=0; i<lenght-j-1; i++)
        {
            if (a[i] <= a[i+1])
            {
                double b = a[i];
                int b_ = index[i];
                a[i] = a[i+1];
                index[i] = index[i+1];

                a[i+1]=b;
                index[i+1] =b_;
            }
        }
    }
}

void bubble_sort_indecies(int *a, int *index, int lenght)
{
    for (int j=0; j<lenght-1; j++)
    {
        for (int i=0; i<lenght-j-1; i++)
        {
            if (a[i] <= a[i+1])
            {
                double b = a[i];
                int b_ = index[i];
                a[i] = a[i+1];
                index[i] = index[i+1];

                a[i+1]=b;
                index[i+1] =b_;
            }
        }
    }
}

void bubble_sort_indecies(double **a, int **index, int lenght, int island)
{
    for (int j=0; j<lenght-1; j++)
    {
        for (int i=0; i<lenght-j-1; i++)
        {
            if (a[island][i] <= a[island][i+1])
            {
                double b = a[island][i];
                int b_ = index[island][i];
                a[island][i] = a[island][i+1];
                index[island][i] = index[island][i+1];

                a[island][i+1]=b;
                index[island][i+1] =b_;
            }
        }
    }
}

void initializePopulation(double **x, double **y, int pop_size, int N, int a, int b)
{
    double ba = b-a;
    for (int i=0; i!=pop_size; i++)
    {
        for (int j=0; j!=N; j++)
        {
            x[i][j] = y[i][j] = RANDOM()*ba+a;
        }
    }
}

void initializeHistory(double ***history_F, double ***history_CR, int H, int M, int island)
{
    for (int i=0; i!=M; i++)
    {
        for (int j=0; j!=H; j++)
        {
            history_F[island][i][j]  = 0.5;
            history_CR[island][i][j] = 0.5;
        }
    }
}

void chooseCrossoverIndecies(int &r1, int &r2, int pbest, int pop_size, int **A, int p, int island)
{
    r1 = RANDOM()*(pop_size-1);
    r2 = RANDOM()*(pop_size-1+A[island][p]);

    while (r1 == r2 || r1 == pbest || r2 == pbest)
    {
        r1 = RANDOM()*(pop_size-1);
        r2 = RANDOM()*(pop_size-1+A[island][p]);
    }
}

void chooseCrossoverIndecies_tournament(int &r1, int &r2, int pbest, int pop_size, int **A, int p, int island, double ***fitness)
{
    r1 = r2 = -1;


    while (r1 == r2 || r1 == pbest || r2 == pbest)
    {
        int t1, t2;
        double fit1, fit2;
        t1=t2=0;
        while (t1==t2)
        {
            t1 = RANDOM()*(pop_size-1);
            t2 = RANDOM()*(pop_size-1);
        }
        fit1=fitness[island][p][t1];
        fit2=fitness[island][p][t2];

        if (fit1<fit2)
        {
            r1 = t1;
        }
        if (fit1>=fit2)
        {
            r1=t2;
        }
        r2 = RANDOM()*(pop_size-1+A[island][p]);
    }
}




void findBestIndex(double ***fitness, int &pbest, int pop_size, int piece_int, int p, int island)
{
    double *f_sort = new double [pop_size];
    int *index_sort = new int [pop_size];

    for (int i=0; i!=pop_size; i++)
    {
        f_sort[i] = fitness[island][p][i];
    }

    int amount = pop_size;
    bubble_sort(f_sort, amount);

    for (int i=0; i!=pop_size; i++)
    {
        for (int j=0; j!=pop_size; j++)
        {
            if (f_sort[i]==fitness[island][p][j])
            {
                index_sort[i]=j;
            }
        }
    }

    pbest = index_sort[int(RANDOM()*piece_int)];

    delete [] f_sort;
    delete [] index_sort;
}

double findMinElement(double *fitness, int pop_size)
{
    int index = 0;
    double minn = fitness[index];

    for (int i=1; i!=pop_size; i++)
    {
        if (fitness[i]<minn)
        {
            minn = fitness[i];
            index = i;
        }
    }

    return minn;
}



void generation_CR(double &CR, double ***H, int r, int p, int island)
{
    CR = randn(H[island][p][r],0.1);
    while (CR>1)
    {
        CR=randn(H[island][p][r],0.1);
    }
    if (CR<0)
    {
        CR=0;
    }
}

void generation_F(double &F, double ***H, int r, int p, int island)
{
    F=randc(H[island][p][r],0.1);
    while (F<-0 || F>1)
    {
        F=randc(H[island][p][r],0.1);
    }
}

void reset_k(int *k, int H, int p)
{
    if (k[p]>=H)
    {
        k[p]=0;
    }
}

void Algorithm_1(double *delta_f, double *w, double *S_CR, double *S_F, double ***HISTORY_CR, double ***HISTORY_F, int *k, int &success, int H, int p, int island)
{

    if (success == 0)
    {
        k[p]++;
        reset_k(k,H,p);
        return;
    }

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int i=0; i!=success; i++)
    {
        sum1+=delta_f[i];
    }

    if (sum1 == 0)
    {
        success = 0;
        k[p]++;
        reset_k(k,H,p);
        return;
    }

    sum1 = 0.0;
    sum2 = 0.0;
    for (int i=0; i!=success; i++)
    {
        sum1+=S_CR[i];
        sum2+=S_F[i];
    }

    if (sum1==0.0 && sum2==0.0)
    {
        k[p]++;
        success = 0;
        reset_k(k,H,p);
        return;
    }

    sum1 = 0.0;
    sum2 = 0.0;
    for (int i=0; i!=success; i++)
    {
        sum1+=S_CR[i];
        sum2+=S_F[i];
    }

    if (sum1==0.0 && sum2!=0.0)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        for (int i=0; i!=success; i++)
        {
            sum1+=delta_f[i];
        }

        for (int i=0; i!=success; i++)
        {
            w[i]=delta_f[i]/sum1;
        }
        sum1 = 0.0;
        sum2 = 0.0;

        for (int i=0; i!=success; i++)
        {
            sum1+=w[i]*S_F[i]*S_F[i];
            sum2+=w[i]*S_F[i];
        }
        double meanF = sum1/sum2;

        HISTORY_F[island][p][k[p]] = meanF;
        if (HISTORY_F[island][p][k[p]] != HISTORY_F[island][p][k[p]])
        {
            HISTORY_F[island][p][k[p]] = 0.5;
        }
        HISTORY_CR[island][p][k[p]] = 0.0;
        k[p]++;
        success=0;
        reset_k(k, H, p);
        return;
    }

    sum1 = 0.0;
    sum2 = 0.0;
    for (int i=0; i!=success; i++)
    {
        sum1+=S_CR[i];
        sum2+=S_F[i];
    }

    if (sum1!=0.0 && sum2==0.0)
    {

        sum2 = 0.0;
        for (int i=0; i!=success; i++)
        {
            sum2+=delta_f[i];
        }
        for (int i=0; i!=success; i++)
        {
            w[i]=delta_f[i]/sum2;
        }
        sum1 = 0.0;
        sum2 = 0.0;

        for (int i=0; i!=success; i++)
        {
            sum1+=w[i]*S_CR[i]*S_CR[i];
            sum2+=w[i]*S_CR[i];
        }
        double meanCR = sum1/sum2;

        HISTORY_CR[island][p][k[p]] = meanCR;
        if (HISTORY_CR[island][p][k[p]] != HISTORY_CR[island][p][k[p]])
        {
            HISTORY_CR[island][p][k[p]]=0.5;
        }
        HISTORY_F[island][p][k[p]] = 0;
        k++;
        success=0;
        reset_k(k,H,p);
        return;
    }


    sum1 = 0.0;
    sum2 = 0.0;
    double sum3 = 0.0;

    for (int i=0; i!=success; i++)
    {
        sum1+=S_CR[i];
        sum2+=S_F[i];
    }

    for (int i=0; i!=success; i++)
    {
        sum3+=delta_f[i];
    }

    if (sum1!=0.0 && sum2!=0.0 && sum3!=0.0 && success!=0)
    {
        sum1 = 0.0;

        for (int i=0; i!=success; i++)
        {
            sum1+=delta_f[i];
        }

        for (int i=0; i!=success; i++)
        {
            w[i]=delta_f[i]/sum1;
        }
        sum1 = 0.0;
        sum2 = 0.0;

        for (int i=0; i!=success; i++)
        {
            sum1+=w[i]*S_F[i]*S_F[i];
            sum2+=w[i]*S_F[i];
        }
        double meanF = sum1/sum2;

        HISTORY_F[island][p][k[p]] = meanF;
        if (HISTORY_F[island][p][k[p]]!=HISTORY_F[island][p][k[p]])
        {
            HISTORY_F[island][p][k[p]] =0.5;
        }

        sum1 = 0.0;
        sum2 = 0.0;

        for (int i=0; i!=success; i++)
        {
            sum1+=w[i]*S_CR[i]*S_CR[i];
            sum2+=w[i]*S_CR[i];
        }

        double meanCR = sum1/sum2;

        HISTORY_CR[island][p][k[p]] = meanCR;
        if (HISTORY_CR[island][p][k[p]] != HISTORY_CR[island][p][k[p]])
        {
            HISTORY_CR[island][p][k[p]] =0.5;
        }
        k[p]++;
        success=0;
        reset_k(k,H,p);
        return;
    }
}


void check_out_borders(double **u, double **population, int i, int N, int a, int b, int *range, int p, int **indeces, int island)
{
    for (int j=range[p]; j!=range[p+1]; j++)
    {
        while (u[i][indeces[island][j]]<a)
        {
            u[i][indeces[island][j]] = a+(abs(a-u[i][indeces[island][j]])/2);
        }
        while (u[i][indeces[island][j]]>b)
        {
            u[i][indeces[island][j]] = b-(abs(b-u[i][indeces[island][j]])/2);
        }
    }
}

void updateArchive(double **archive, double **population, int i, int archive_size, int **A, int *range, int p, int **indeces, int island)
{
    if (A[island][p]>=archive_size)
    {
        int index = RANDOM()*(archive_size-1);
        for (int j=range[p]; j!=range[p+1]; j++)
        {
            archive[index][indeces[island][j]] = population[i][indeces[island][j]];
        }
        return;
    }

    if (A[island][p]<archive_size)
    {
        for (int j=range[p]; j!=range[p+1]; j++)
        {
            archive[A[island][p]][indeces[island][j]] = population[i][indeces[island][j]];
        }
        A[island][p]++;

        return;
    }
}


double mean_stat(double **x, int gener, int R)
{

    double sum;
    sum = 0.0;

    for (int i = 0; i<R; i++)
    {
        sum += x[i][gener];
    }

    sum = sum / (double)R;

    return sum;
}


double min_stat(double **x, int gener, int R)
{

    double sum;

    sum = x[0][gener];

    for (int i = 1; i<R; i++)
    {
        if (x[i][gener]<sum)
        {
            sum = x[i][gener];
        }
    }
    return sum;
}

double max_stat(double **x, int gener, int R)
{

    double sum;

    sum = x[0][gener];

    for (int i = 1; i<R; i++)
    {
        if (x[i][gener]>sum)
        {
            sum = x[i][gener];
        }
    }


    return sum;
}

double median_stat(double **x, int gener, int R)
{
    double answ;
    answ = 0.0;
    double *fitness = new double[R];

    for (int i = 0; i<R; i++)
    {
        fitness[i] = x[i][gener];
    }

    quickSort(fitness, 0, R-1);

    R = (R - 1) / 2;
    answ = fitness[R];

    delete[]fitness;

    return answ;
}

double stddev_stat(double **x, int gener, int R, double M)
{
    double sum;
    sum = 0.0;

    for (int i = 0; i<R; i++)
    {
        sum += pow((M - x[i][gener]), 2);
    }

    R = R - 1;
    sum = sum / (double)R;
    sum = sqrt(sum);

    return sum;

}


void rnd_indecies(int **vec, int pop_size, int M, int island)
{
    for (int i=0; i!=island; i++)
    {
        for (int j=0; j!=M; j++)
        {
            vec[i][j] = RANDOM()*(pop_size-1);
        }
    }
}

void indecesSuccession(int *x, int N)
{
    for (int i = 0; i != N; i++)
    {
        x[i] = i;
    }
}

void indecesSuccession(int **x, int pool_size, int N)
{
    for (int i=0; i!=pool_size; i++)
    {
        for (int j=0; j!=N; j++)
        {
            x[i][j] = j;
        }
    }
}


void randperm(int *x, int N)
{
    int a, temp;
    for (int i = 0; i != N; i++)
    {
        temp = x[i];
        a = RANDOM()*(N - 1);
        swap(x[i], x[a]);
    }
}
void randperm(int **x, int N, int index)
{
    int a, temp;
    for (int i = 0; i != N; i++)
    {
        temp = x[index][i];
        a = RANDOM()*(N - 1);
        swap(x[index][i], x[index][a]);
    }
}

void randperm_in_all_islands(int **x, int N, int pool_size)
{
    int a, temp;
    temp =0;
    int index =0;
    for (int i = 0; i != N; i++)
    {
        temp = x[index][i];
        a = RANDOM()*(N - 1);
        swap(x[index][i], x[index][a]);
    }
    for (int i=index+1;i!=pool_size;i++)
    {
        for (int j=0;j!=N;j++)
        {
            x[i][j] = x[0][j];
        }
    }
}

void randperm_island(int **indices, int N, int pool_size)
{
    int a, temp;

    for (int i = 0; i != N; i++)
    {
        temp = indices[0][i];
        a = RANDOM()*(N - 1);
        swap(indices[0][i], indices[0][a]);
    }
    for (int k=1; k!=pool_size; k++)
    {
        for (int i = 0; i != N; i++)
        {
            indices[k][i] = indices[0][i];
        }
    }
}
void randperm_island_for_all(int **indices, int N, int pool_size)
{
    int a, temp;

    for (int i = 0; i != N; i++)
    {
        temp = indices[0][i];
        a = RANDOM()*(N - 1);
        swap(indices[0][i], indices[0][a]);
    }
    for (int k=1; k!=pool_size; k++)
    {
        for (int i = 0; i != N; i++)
        {
            indices[k][i] = indices[0][i];
        }
    }
}
void find_best_part_index (int **best_indecies, double ***cc_fitness, int p, int pop_size, int island)
{
    best_indecies[island][p]=0;
    double min_fitness = cc_fitness[island][p][0];
    for (int j=0; j!=pop_size; j++)
    {
        if (cc_fitness[island][p][j]<min_fitness)
        {
            min_fitness = cc_fitness[island][p][j];
            best_indecies[island][p] = j;
        }
    }
}

double find_best_fitness_value (double ***fitness_cc, int M, int *pop_size, int island)
{
    double minn = fitness_cc[island][0][0];

    for (int i=0; i!=M; i++)
    {
        for (int j=0; j!=pop_size[island]; j++)
        {
            if (fitness_cc[island][i][j]<minn)
            {
                minn = fitness_cc[island][i][j];
            }
        }
    }
    return minn;
}


void filling_zeros(double *arr, int N)
{
    for (int i=0; i!=N; i++)
    {
        arr[i]=0.0;
    }
}

void filling_zeros_A(int **arr, int N, int pool_size)
{
    for (int i=0; i!=N; i++)
    {
        for (int j=0; j!=pool_size; j++)
            arr[j][i]=0.0;
    }
}
void filling_zeros_A(double **arr, int N, int pool_size)
{
    for (int i=0; i!=N; i++)
    {
        for (int j=0; j!=pool_size; j++)
            arr[j][i]=0.0;
    }
}

void filling_zeros_A(double ***arr, int a, int b, int c)
{
    for (int i=0; i!=a; i++)
    {
        for (int j=0; j!=b; j++)
        {
            for (int k=0; k!=c; k++)
            {
                arr[i][j][k]=0.0;
            }
        }
    }
}


void filling_zeros(double**arr, int M, int N)
{
    for (int i=0; i!=M; i++)
    {
        for (int j=0; j!=N; j++)
        {
            arr[i][j]=0.0;
        }

    }
}

void filling_number_double_arr(double **arr, int N, int M, double number)
{
    for (int i=0; i!=N; i++)
    {
        for (int j=0; j!= M; j++)
            arr[i][j]=number;
    }
}

void filling_zeros(int *arr, int N)
{
    for (int i=0; i!=N; i++)
    {
        arr[i]=0.0;
    }
}

void building_solution(int *range, double *solution, int **indeces, int **cc_best_individual_index, double **population, int i, int p, int M, int island)
{
    for (int j=range[p]; j!=range[p+1]; j++)
    {
        solution[indeces[island][j]] = population[i][indeces[island][j]];
    }
    for (int p_cc = 0; p_cc < M; p_cc++)
    {
        if (p != p_cc)
        {
            for (int j = range[p_cc]; j < range[p_cc + 1]; j++)
            {
                solution[indeces[island][j]] = population[cc_best_individual_index[island][p_cc]][indeces[island][j]];
            }
        }
    }
}

void building_solution_U(int *range, double *solution, double **u, int **indeces, int **cc_best_individual_index, double **population, int i, int p, int M, int island)
{
    for (int j=range[p]; j!=range[p+1]; j++)
    {
        solution[indeces[island][j]] = u[i][indeces[island][j]];
    }

    for (int p_cc = 0; p_cc < M; p_cc++)
    {
        if (p != p_cc)
        {
            for (int j = range[p_cc]; j < range[p_cc + 1]; j++)
            {
                solution[indeces[island][j]] = population[cc_best_individual_index[island][p_cc]][indeces[island][j]];
            }
        }
    }
}

void copying_from_first_to_second (double *arr_1, double *arr_2, int N)
{
    for (int i=0; i!=N; i++)
    {
        arr_2[i] = arr_1[i];
    }
}

void copying_from_first_to_second_2D (double **arr_1, double **arr_2, int M, int N)
{
    for (int i=0; i!=M; i++)
    {
        for (int j=0; j!=N; j++)
        {
            arr_2[i][j] = arr_1[i][j];
        }
    }

}

void information_output (string name_of_func, double THE_GLOBAL_BEST_solution,double THE_BEST_solution, int FEV, int z, int R, int M, int pop_size, int archive_size, int H, double *** HISTORY_CR, double ***HISTORY_F, int island, int restarts)
{
    std::cout<<name_of_func<<endl;
    std::cout<<"Global best solution value: "<<THE_GLOBAL_BEST_solution<<endl;
    std::cout<<"Current best solution value: "<<THE_BEST_solution<<endl;
    std::cout<<"Current FEV: "<<FEV<<endl;
    std::cout<<z+1<<" out of "<<R<<" RUNS"<<endl;
    std::cout<<"Number of subcomponents: "<<M<<endl;
    std::cout<<"Current pop_size: "<<pop_size<<endl;
    std::cout<<"The numner of Restarts: "<<restarts<<endl;
    cout << "==========================================================" << endl<<endl;

}

void GENERATING_TRIAL_VECTOR (int i, double ***fitness_cc, int &pbest, int pop_size, double piece, int p, int *r, int H, double *F, double *CR, double ***HISTORY_CR,
                              double ***HISTORY_F, int r1, int r2, int **A, int *range, double **u, double **population, int **indeces, double **archive, int island, bool tournament_selection)
{
    int piece_int = pop_size*piece;
    findBestIndex(fitness_cc, pbest, pop_size, piece_int, p, island);
    r[i] = RANDOM() * (H-1);
    generation_CR(CR[i], HISTORY_CR, r[i], p, island);
    generation_F(F[i], HISTORY_F, r[i], p, island);
    if (tournament_selection == 0)
    {
        chooseCrossoverIndecies(r1, r2, pbest, pop_size, A, p, island);
    }
    if (tournament_selection == 1)
    {
        chooseCrossoverIndecies_tournament(r1, r2, pbest, pop_size, A, p, island, fitness_cc);
    }

    if (r2<pop_size)
    {
        for (int j=range[p]; j!=range[p+1]; j++)
        {
            u[i][indeces[island][j]] = population[i][indeces[island][j]]+F[i]*(population[pbest][indeces[island][j]]-population[i][indeces[island][j]])+F[i]*(population[r1][indeces[island][j]]-population[r2][indeces[island][j]]);
        }
    }

    if (r2>=pop_size)
    {
        r2-=pop_size;
        for (int j=range[p]; j!=range[p+1]; j++)
        {
            u[i][indeces[island][j]] = population[i][indeces[island][j]]+F[i]*(population[pbest][indeces[island][j]]-population[i][indeces[island][j]])+F[i]*(population[r1][indeces[island][j]]-archive[r2][indeces[island][j]]);
        }
    }
    int jrand = RANDOM()*((range[p+1]-range[p])+range[p]);

    for (int j=range[p]; j!=range[p+1]; j++)
    {
        if(RANDOM()<=CR[i] || j==jrand)
        {
        }
        else
        {
            u[i][indeces[island][j]] = population[i][indeces[island][j]];
        }
    }
}

void UPDATING_INFO(double **archive, double **population, int i, int archive_size, int **A, int *range, int p, int **indeces,
                   double ***fitness_cc, double ***fitness_cc_new, double test_function, double *delta_f, double *S_F, double *F, double *S_CR, double *CR, int &success, int island)
{
    updateArchive(archive, population, i, archive_size, A, range, p, indeces, island);
    fitness_cc_new[island][p][i] = test_function;
    delta_f[success] = sqrt ((test_function-fitness_cc[island][p][i])*(test_function-fitness_cc[island][p][i]));
    S_F[success] = F[i];
    S_CR[success] = CR[i];
    success++;
}


void max_CC_parts (int *parts, int &max_part, int amount)
{
    max_part = 0;
    for (int i=0; i!=amount; i++)
    {
        if (max_part<parts[i])
        {
            max_part = parts[i];
        }
    }
}

int finding_the_best_index(int pool_size, double *best_pool_fitness)
{
    int the_best_index = 0;
    double min_value = best_pool_fitness[0];

    for (int i=0; i!=pool_size; i++)
    {
        if (min_value>best_pool_fitness[i])
        {
            min_value = best_pool_fitness[i];
            the_best_index = i;
        }
    }

    return the_best_index;
}

void recalculating_range_interval(int *range, int M, int N, int S)
{
    range[0] = 0;
    range[M] = N;
    S = N / M;
    for (int i = 1; i < M; i++)
    {
        range[i] = range[i - 1] + S;
    }
}


void migration(double ***population, double ***population_new, double *best_island_fitness, double **best_island_solution, int *pop_size, int pool_size, int *island_pool, int *indeces, int N, int **cc_best_individual_index, int &best_index, double ***fitness_cc)
{
    double best_solution = 1e300;

    for (int i=0; i!=pool_size; i++)
    {
        if (best_solution>best_island_fitness[i])
        {
            best_solution = best_island_fitness[i];
            best_index = i;
        }
    }

    for (int i=0; i!=pool_size; i++)
    {
        int rnd_index = RANDOM()* (pop_size[i]-1);

        if (i != best_index)
        {
            for (int j=0; j!=N; j++)
            {
                population[i][rnd_index][indeces[j]] = population_new[i][rnd_index][indeces[j]] = best_island_solution[best_index][indeces[j]];
            }
            for (int k=0; k!=island_pool[i]; k++)
            {
                cc_best_individual_index[i][k] = rnd_index;
            }

            for(int k=0; k!=island_pool[i]; k++)
            {
                fitness_cc[i][k][rnd_index] = best_solution;

            }

        }
    }
}

double mean_x(double **population, int pop_size, int N, int index)
{
    double sum = 0.0;
    for (int i=0; i!=pop_size; i++)
    {
        sum+=population[i][index];
    }
    return sum/pop_size;
}

double DI_calculation (double **population, int pop_size, int N)
{
    double DI = 0.0;

    for (int i=0; i!=pop_size; i++)
    {
        for(int j=0; j!=N; j++)
        {
            double x_aver = mean_x(population, pop_size, N, j);
            DI+=(population[i][j]-x_aver)*(population[i][j]-x_aver);
        }
    }

    return sqrt(DI/pop_size);
}

void increase_population (double ***population, int *pop_size,  int island, int N, double a, double b)

{
    pop_size[island]++;
    for (int j=0; j!=N; j++)
    {
        population[island][pop_size[island]-1][j] = RANDOM()*(b-a)+a;
    }


}


void borbers_check(double *newsol, int N, int a, int b)
{
    for (int i=0; i!=N; i++)
    {
        while (newsol[i]<a)
        {
            newsol[i]=a+(abs(newsol[i]-a))/2;
        }
        while (newsol[i]>b)
        {
            newsol[i]=b-(abs(newsol[i]-b))/2;
        }

    }
}

void mix_sequence(int *arr, int lenght)
{
    for (int i=0; i!=lenght; i++)
    {
        int a = RANDOM()*lenght;
        swap(arr[a],arr[i]);
    }


}


void DI_Algorithm (double **population, double **population_new, int &pop_size, int N, int DI_init, int FEV_global, int FEV, int pop_size_min,
                   int pop_size_max, int a, int b, int *island_pool, int *range, int S, double ***fitness_cc,
                   double ***fitness_cc_new, int island)
{

    double DI = DI_calculation(population, pop_size, N);
    double RD = DI/DI_init;
    double RFES = (double)(FEV_global-FEV)/(double)FEV_global;
    double rRD = 1-RFES/0.9;

    if (rRD<0.0)
    {
        rRD = 0.0;
    }

    int new_pop_size_increased = pop_size+1;

    if (RD<0.9*rRD && new_pop_size_increased<=pop_size_max)
    {
        for (int i=0; i!=N; i++)
        {
            population[new_pop_size_increased-1][i] = population_new[new_pop_size_increased-1][i] = RANDOM()*(b-a)+a;
        }
        int M = island_pool[island];
        for (int p=0;p!=M;p++)
        fitness_cc[island][p][new_pop_size_increased-1] = fitness_cc_new[island][p][new_pop_size_increased-1] = 1e300;

        pop_size=new_pop_size_increased;
    }

    int new_pop_size_downcreased = pop_size-1;

    if ((RD>1.1*rRD) && new_pop_size_downcreased  >= pop_size_min)
    {
        int M = island_pool[island];
        recalculating_range_interval(range, M, N, S);

        for (int p=0; p!=M; p++)
        {
            int worst_index = 0;
            double worst_fitness = fitness_cc[island][p][0];
            for (int i=0; i!=pop_size; i++)
            {
                if (fitness_cc[island][p][i]>worst_fitness)
                {
                    worst_index = i;
                    worst_fitness = fitness_cc[island][p][i];
                }
            }
            for (int i=worst_index; i!=pop_size-1; i++)
            {
                for (int k=range[p]; k!=range[p+1]; k++)
                {
                    population[i][k] = population_new[i][k] = population[i+1][k];
                }
                fitness_cc[island][p][i] = fitness_cc_new[island][p][i] = fitness_cc[island][p][i+1];
            }
        }
        pop_size=new_pop_size_downcreased;
    }
}

void LS_SR_update(int i, int *LS1_indecies, double *LS1_improve, double *SR, int a, int b)
{
    if (LS1_improve[LS1_indecies[i]] == 0.0)
    {
        SR[LS1_indecies[i]]=SR[LS1_indecies[i]]*0.5;
    }

    if (SR[LS1_indecies[i]]<1e-18)
    {
        SR[LS1_indecies[i]] = 0.4*(b-a);
    }
}

void Improvement (double improvement, double THE_BEST_solution, double threshold, int &num_worse, bool &trigger_restart)
{

    double sol;
    if (improvement == 0.0)
    {
        sol = 1.0;
    }
    else
    {
        sol = (improvement-THE_BEST_solution)/improvement;

    }


    if (sol>threshold)
    {
        num_worse = 0;
    }
    if (sol<=threshold)
    {
        num_worse++;
    }

    if (num_worse >=3)
    {
        trigger_restart = 1;
        num_worse = 0;
    }
}

double median_island (double threshold, int pop_size, int island, int M, double ***fitness_cc)
{

    double *arr = new double [pop_size*M];
    int counter = 0;
    for (int i=0; i!=M; i++)
    {
        for (int j=0; j!=pop_size; j++)
        {
            arr[counter] = fitness_cc[island][i][j];
            counter++;
        }
    }
    bubble_sort(arr, counter);

    int index = threshold * pop_size/2;
    double otv = arr[index];
    delete []arr;
    return otv;
}

double median_all_islands (double threshold, int pop_size, double ***fitness_cc, int pool_size, int *island_pool)
{
    int total_size_arr = 0;
    for (int i=0;i!=pool_size;i++)
    {
        total_size_arr+=pop_size*island_pool[i];
    }

    double *arr = new double [total_size_arr];
    int counter = 0;
    for (int island=0; island!=pool_size; island++)
    {
        for (int i=0; i!=island_pool[island]; i++)
        {
            for (int j=0; j!=pop_size; j++)
            {
                arr[counter] = fitness_cc[island][i][j];
                counter++;
            }
        }
    }
    bubble_sort(arr, counter);

    int index = threshold * pop_size/2;
    double otv = arr[index];
    delete []arr;
    return otv;
}

void FEV_island_calculating(int *FEV_island_generations, int *FEV_island_generations_previous, int min_generations, int max_generations, double *island_performance, int pool_size)
{

    int best_index = 0;
    double best_performance = island_performance[0];
    for (int i=0;i<=pool_size;i++)
    {
        if (island_performance[i]>best_performance)
        {
            best_index = i;
            best_performance = island_performance[i];
        }
    }
    for (int i=0;i<=pool_size;i++)
    {
        FEV_island_generations[i] = FEV_island_generations_previous[i];
    }

    int generation_pool = 0;
    for (int i=0;i!=pool_size;i++)
    {
        if (FEV_island_generations[i]-1 >=min_generations)
        {
            FEV_island_generations[i]--;
            generation_pool++;
        }
    }

    FEV_island_generations[best_index]+=generation_pool;

    for (int i=0;i!=pool_size;i++)
    {
        FEV_island_generations_previous[i] = FEV_island_generations[i];
    }

}
