#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(void)
{
	static unsigned int seed = 0;
    ++seed;

    srand((unsigned)time(NULL)+seed*seed);

    printf("%f\n", rand()*1.0/(RAND_MAX)); /*TODO: 检查对不对??*/
}
