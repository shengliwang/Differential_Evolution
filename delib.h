#ifndef __DELIB_H__
#define __DELIB_H__


typedef double (*fitnessFunc_t)(double * argv, int argc);

typedef struct DE_args
{
    fitnessFunc_t fitnessFunc;
    double F;   /* 变异幅度 */
    double CR;  /* 交叉率 */
    unsigned int NP; /* number of Population, 种群大小*/
    unsigned int ND; /* nubmer of Dimensions, 个体的大小,即解空间维数*/
}DE_INIT_ARG;

int delib_gen_one_step(double * result);
int delib_init(DE_INIT_ARG * arg);
int delib_deinit(void);





#endif
