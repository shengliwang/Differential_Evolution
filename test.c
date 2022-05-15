#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "delib.h"

#define PI (acos(-1))


/*Bent Cigar Function: f=x1*x1 + 10^6(Σxi2)*/
double Bent_Cigar_Func(double * argv, int argc)
{
    if (argc < 2) /*此函数至少要两个变量*/
    {
        printf("Bent Cigar Function < 2\n");
        return 0;
    }

    double tmp = 0;
    for (int i = 1; i < argc; ++i)
    {
        tmp += argv[i] * argv[i];
    }

    return argv[0]*argv[0] + 1000000*tmp;
}

/*抛物线函数*/
double paowuxian_func(double * args, int argc)
{
    double fitness = 0;
    for (int i = 0; i < argc; ++i)
    {
        fitness += (args[i]+1) * (args[i]+1);
    }

    return fitness+2;
}

/* Rastrigin's Func*/
double Rastrigin_Func(double * argv, int argc)
{
    double tmp = 0;
    for (int i = 0; i < argc; ++i)
    {
        tmp += (argv[i]*argv[i] - 10*cos(2*PI*argv[i]) + 10);
    }

    return tmp;
}

/*Rosenbrock's Function*/
double Rosenbrock_Function(double * x, int argc)
{
    double tmp = 0;
    if (argc < 2)
    {
        printf("warning: (argc < 2)\n");
        return 0;
    }

    for (int i = 0; i < argc -1; ++i)
    {
        tmp += (100*(x[i]*x[i]-x[i+1])*(x[i]*x[i]-x[i+1])
                + (x[i]-1)*(x[i]-1));
    }
    return tmp;
}


int main(void)
{
    /*测试1: 测试Bent Cigar Function*/
    DE_INIT_ARG de_args;
    de_args.fitnessFunc = Rastrigin_Func;
    de_args.CR = 0.9;
    de_args.F = 0.5;
    de_args.NP = 1000;
    de_args.ND = 2;

    double * result = malloc(sizeof(double)*de_args.ND);

    
    delib_init(&de_args);

    int step = 0;
    while(step <= 1000)
    {
        delib_gen_one_step(result);

        printf("step %d: x=[", step);
        for (int i = 0; i < de_args.ND; ++i)
        {
            printf("%f ", result[i]);
        }
        printf("]. result=%f\n", de_args.fitnessFunc(result, de_args.ND));
        ++step;
    }



    delib_deinit();
}
