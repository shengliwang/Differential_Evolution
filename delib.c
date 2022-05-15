/*
 * 本文件实现DE算法(Differential Evolution)的算法库(C语言版本)
 * 作者: 王胜利
 * Email: vicking@mail.ustc.edu.cn
 * 时间: 2022-5-19
 * 参考: https://github.com/zhangzhiyong1989/Differential-Evolution/blob/master/DE.py
 * 说明: 变量前的前缀表示当前变量的变量类型.
 *       如uiPopulationN表示这个变量为unsigned int类型
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <time.h>
 #include <string.h>

/*生成的随机整数n: min <= n <= max*/
static int myrand_int(int min, int max)
{
    if (min >= max)
    {
        printf("warning: min >= max\n");
    }

    static unsigned int seed = 0;
    ++seed;

    srand((unsigned)time(NULL)+seed*seed);
    int n = rand();

    return min + n%(max - min + 1);
}

/*产生0-1之间的随机数*/
static double myrand_double()
{
    static unsigned int seed = 0;
    ++seed;

    srand((unsigned)time(NULL)+seed*seed);

    return rand()*1.0/(RAND_MAX);
}

static double ** malloc_matrix(unsigned int uiRows, unsigned int uiCols)
{
    int loop;
    double ** tmp = malloc(sizeof(double *)*uiRows);
    if (NULL == tmp)
    {
        return NULL;
    }

    for (loop = 0; loop < uiRows; ++loop)
    {
        tmp[loop] = malloc(sizeof(double) * uiCols);
        if (NULL == tmp[loop])
        {
            return NULL;
        }
        memset(tmp[loop], 0, sizeof(double) * uiCols);
    }

    return tmp;
}

static void free_matrix(double ** mat, unsigned int uiRows)
{
    int loop = 0;
    for (; loop < uiRows; ++loop)
    {
        free(mat[loop]);
    }

    free(mat);
}


double calFitness(const double * args, unsigned int argc)
{
    double fitness = 0;
    for (int i = 0; i < argc; ++i)
    {
        fitness += (args[i]+1) * (args[i]+1);
    }

    return fitness+2;
}


/*
 * 执行变异的函数
 * [入参]matrixPopulation: 表示种群的矩阵
 * [入参]uiPopulationN: 种群中个体的数量
 * [入参]uiIndividualSize: 每个个体的大小(即解空间的维数)
 * [入参]dF: 变异幅度
 * [返回值]: 返回 根据种群变异的临时种群
 * Note: 返回值为malloc申请获得,注意要释放内存.
*/
double ** mutation(double ** matrixPopulation, 
        unsigned int uiPopulationN, unsigned int uiIndividualSize, double dF)
{
    if (NULL == matrixPopulation)
    {
        printf("matrixPopulation is NULL\n");
        return NULL;
    }

    double ** matrixPopulationTmp = malloc_matrix(uiPopulationN, uiIndividualSize);
    if (NULL == matrixPopulationTmp)
    {
        printf("matrixPopulationTmp is NULL\n");
        return NULL;
    }

    /*产生新的变异个体*/
    for (int i = 0; i < uiPopulationN; ++i)
    {
        int r1 = 0, r2 = 0, r3 = 0;
        while(r1 == i || r2 == i || r3 == i || r1 == r2 || r1 == r3 || r2 == r3)
        {
            r1 = myrand_int(0, uiPopulationN-1);
            r2 = myrand_int(0, uiPopulationN-1);
            r3 = myrand_int(0, uiPopulationN-1);
        }

        for (int j = 0; j < uiIndividualSize; ++j)
        {
            matrixPopulationTmp[i][j] =
              matrixPopulation[r1][j] + dF * (matrixPopulation[r2][j] - matrixPopulation[r3][j]);
        }
    }
    
    return matrixPopulationTmp;
}


/*
 * 执行交叉的函数
 * [入参]matrixPopulation: 表示种群的矩阵
 * [入参]matrixPopulationTmp:调用mutation函数产生的临时种群
 * [入参]uiPopulationN: 种群中个体的数量
 * [入参]uiIndividualSize: 每个个体的大小(即解空间的维数)
 * [入参]dCR:交叉率
 * [返回值]: 返回 种群和mutation产生的临时种群 交叉后的种群
 * Note: 返回值为malloc申请获得,注意要释放内存.
*/
double ** crossover(double ** matrixPopulation, double ** matrixPopulationTmp,
        unsigned int uiPopulationN, unsigned int uiIndividualSize, double dCR)
{
    if (NULL == matrixPopulation || NULL == matrixPopulationTmp)
    {
        printf("(NULL == matrixPopulation || NULL == matrixPopulationTmp)\n");
        return NULL;
    }

    double ** matrixPopulationCrossOver = malloc_matrix(uiPopulationN, uiIndividualSize);
    if (NULL == matrixPopulationCrossOver)
    {
        printf("matrixPopulationCrossOver is NULL\n");
        return NULL;
    }

    for (int i = 0; i < uiPopulationN; ++ i)
    {
        for (int j = 0; j < uiIndividualSize; ++j)
        {
            double r = myrand_double();
            if (r <= dCR)
            {
                matrixPopulationCrossOver[i][j] = matrixPopulationTmp[i][j];
            }
            else
            {
                matrixPopulationCrossOver[i][j] = matrixPopulation[i][j];
            }
        }
    }

    return matrixPopulationCrossOver;
}

/*
 * 选择函数: 选择出下一代的函数.
 * [入出参]matrixPopulation,表示种群的矩阵(此函数会更新此矩阵)
 * [入  参]  matrixPopulationCrossOver,由crossover交叉后产生的临时矩阵
 * [入  参]  uiPopulationN: 种群中个体的数量
 * [入  参]  uiIndividualSize: 每个个体的大小(即解空间的维数)
 * [入出参]fitnessval,适应值向量(此函数会更新此向量)
*/
void selection(double ** matrixPopulation, double ** matrixPopulationCrossOver,
        unsigned int uiPopulationN, unsigned int uiIndividualSize,double *fitnessval)
{
    if (NULL == matrixPopulation || NULL == matrixPopulationCrossOver || NULL == fitnessval)
    {
        printf("(NULL == matrixPopulation || NULL == matrixPopulationCrossOver)");
        return;
    }
    double * fitnessval_crossover = malloc(sizeof(double)*uiPopulationN);

    for (int i = 0; i < uiPopulationN; ++i)
    {
        fitnessval_crossover[i] = calFitness(matrixPopulationCrossOver[i], uiIndividualSize);
        printf("fitnessval_crossover[%d]=%f, fitnessval[%d]=%f\n",
                        i, fitnessval_crossover[i],
                        i, fitnessval[i]);
        if (fitnessval_crossover[i] < fitnessval[i])
        {
            for (int j = 0; j < uiIndividualSize; ++j)
            {
                matrixPopulation[i][j] = matrixPopulationCrossOver[i][j];
            }
            fitnessval[i] = fitnessval_crossover[i];
        }
    }
    free(fitnessval_crossover);
}


void saveBest(double * fitnessVal, unsigned int uiPopulationN)
{
    if (NULL == fitnessVal)
    {
        printf("(NULL == fitnessVal)");
        return ;
    }
    
    int tmp = 0;

    for (int i = 1; i < uiPopulationN; ++i)
    {
        if (fitnessVal[i] < fitnessVal[tmp])
        {
            tmp = i;
        }
    }

    printf("best val: %f\n", fitnessVal[tmp]);
}


int main(void)
{
    int xMin = -10;
    int xMax = 10;
    /*初始化种群, 1000个个体,解空间为10维*/
    unsigned int NP = 10;
    unsigned int ND = 2; /*number of dimensions*/
    double ** matrixPopulation = malloc_matrix(NP, ND);
    for (int i = 0; i < NP; ++i)
    {
        for (int j = 0; j < ND; ++j)
        {
            matrixPopulation[i][j] = xMin + myrand_double() * (xMax - xMin);
        }
    }

    for (int i = 0; i < NP; ++i)
    {
        printf("indival%d: ", i);
        for (int j = 0; j < ND; ++j)
        {
            printf("%f ", matrixPopulation[i][j]);
        }
        printf("\n");
    }

    /*计算适应值*/
    double * fitnessVal = malloc(sizeof(double)*NP);
    for (int i = 0; i < NP; ++i)
    {
        fitnessVal[i] = calFitness(matrixPopulation[i], ND);
    }

    printf("fitness val:");
    for (int i = 0; i < NP; ++i)
    {
        printf("%f ", fitnessVal[i]);
    }
    printf("\n");

    /*开始演化计算*/
    int gen = 0; 
    double F = 0.5;
    double CR = 0.9;
    
    while (gen <= 1000) /*计算1000代*/
    {
        printf("generation %d\t: ", gen); fflush(stdout);
         printf("\nmatrixPopulation:\n");
        for (int i = 0; i < NP; ++i)
    {
        printf("indival%d: ", i);
        for (int j = 0; j < ND; ++j)
        {
            printf("%f ", matrixPopulation[i][j]);
        }
        printf("\n");
    }
        double ** matrixPopulationMutation = mutation(matrixPopulation, NP, ND, F);
        printf("\nmatrixPopulationMutation:\n");
        for (int i = 0; i < NP; ++i)
    {
        printf("indival%d: ", i);
        for (int j = 0; j < ND; ++j)
        {
            printf("%f ", matrixPopulationMutation[i][j]);
        }
        printf("\n");
    }
        double ** matrixPopulationCrossOver = crossover(matrixPopulation, matrixPopulationMutation,
                                                    NP, ND,CR);
        printf("\nmatrixPopulationCrossOver:\n");
        for (int i = 0; i < NP; ++i)
    {
        printf("indival%d: ", i);
        for (int j = 0; j < ND; ++j)
        {
            printf("%f ", matrixPopulationCrossOver[i][j]);
        }
        printf("\n");
    }
        printf("fitness val1:");
    for (int i = 0; i < NP; ++i)
    {
        printf("%f ", fitnessVal[i]);
    }
    printf("\n");
        selection(matrixPopulation, matrixPopulationCrossOver, NP, ND, fitnessVal);
        printf("fitness val2:");
    for (int i = 0; i < NP; ++i)
    {
        printf("%f ", fitnessVal[i]);
    }
    printf("\n");
        free_matrix(matrixPopulationMutation, NP);
        free_matrix(matrixPopulationCrossOver, NP);
        saveBest(fitnessVal, NP);
        ++gen;
    }
}

