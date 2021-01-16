//function for implementing crosscorelation

#include <android/log.h>
#include "GLD.h"

//function for implementing crosscorrelation
double *pmxcorr(double *x, double *y, int maxlag, int frameSize)
{
    int kk = frameSize-maxlag+1;
    double* l;
    double* rxy;
    double* pm_rxy;

    l = (double*)calloc((2 * kk + 1), sizeof(double));

    int p = 0;
    for (int i = kk; i >= 0; i--)
    {
        l[p] = -1 * i;
        //printf("l[%d], %d\n", p, -1 * i);
        p = p + 1;

    }

    for (int ii = 1; ii <= kk; ii++)
    {
        l[p] = ii;
        //printf("l[%d], %d\n", p, ii);
        p = p + 1;
    }

    rxy = (double*)calloc((2 * kk + 1), sizeof(double));

    for (int g = l[0]; g <= l[2 * kk]; g++)
    {
        for (int i = 0; i < kk; i++)
        {
            if (((i + 1) <= 0) || (((i + 1) - g) <= 0) || (((i + 1) - g) > kk))
            {
                rxy[g + kk] = rxy[g + kk];
                //printf("rxy[%d], %.32lf\n", g + kk, rxy[g + kk]);
            }
            else
            {
                rxy[g + kk] = rxy[g + kk] + x[i] * y[i - g];

                //printf("x[%d], %.32lf\n", i, x[i]);
                //printf("y[%d], %.32lf\n", i - g, y[i - g]);

                //printf("rxy[%d], %.32lf\n", g + kk, rxy[g + kk]);
            }
        }
    }

    pm_rxy = (double*)calloc((2 * maxlag + 1), sizeof(double));

    int d = 0;

    for (int u = floor((2 * kk - 1) / 2) - maxlag + 1; u <= floor((2 * kk - 1) / 2) + maxlag + 1; u++)
    {

        pm_rxy[d] = rxy[u];

        //printf("rxy[%d], %.32lf\n", u, rxy[u]);
        //printf("pm_rxy[%d], %.32lf\n", d, pm_rxy[d]);

        d = d + 1;
    }

    return pm_rxy;
}


//function for finding energy
double energy(double* x,  int order)
{
    double sum = 0;

    for (int i = 0; i <= order; i++)
    {
        sum = sum + x[i] * x[i];

    }
    return sum;
}


//function for implementing generalized levinson durbin algorithm
double *GLD(double *mic_in, double *noise, int maxlag,int frameSize)
{
    double* crosscorr;
    double* yn_corr;
    double* E;
    double* K;
    double*a;
    double*b;
    double*pre_a;
    double*pre_b;
    double*rxx;
    double*R;
    double*R_reverse;
    int i, j, M;
    int p = maxlag + 1;
    double E0;
    double temp1;
    double temp2;
    double* imp;


    E = (double*)calloc(maxlag + 1, sizeof(double));
    K = (double*)calloc(maxlag + 1, sizeof(double));
    a = (double*)calloc(maxlag + 1, sizeof(double));
    b = (double*)calloc(maxlag + 1, sizeof(double));
    pre_a = (double*)calloc(maxlag + 1, sizeof(double));
    pre_b = (double*)calloc(maxlag + 1, sizeof(double));
    rxx = (double*)calloc(maxlag + 3, sizeof(double));
    R = (double*)calloc(maxlag + 1, sizeof(double));
    R_reverse = (double*)calloc(maxlag + 1, sizeof(double));
    yn_corr = (double*)calloc((maxlag + 1), sizeof(double));
    imp = (double*)calloc(maxlag + 1, sizeof(double));

    /*
    clock_t begin = clock();

    crosscorr = pmxcorr(mic_in, noise, maxlag , frameSize);

    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    __android_log_print(ANDROID_LOG_INFO, "k", "meeeeeeeeeeeeeeeeeeeeeeerrrrrrrrrrrrrrraaaaaaaaaaaaaaaaaaaaaatiming = %1.7g\n", time_spent);
*/
  /*  clock_t start, end;
    double cpu_time_used;

    start = clock();

    crosscorr = pmxcorr(mic_in, noise, maxlag , frameSize);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    __android_log_print(ANDROID_LOG_INFO, "k", "meeeeeeeeeeeeeeeeeeeeeeerrrrrrrrrrrrrrraaaaaaaaaaaaaaaaaaaaaatimingcrosscorr = %1.7g\n", cpu_time_used);
*/

    crosscorr = pmxcorr(mic_in, noise, maxlag , frameSize);


    int ddd = 0;
    for (j = 0; j < maxlag + 1; j++)
    {
        yn_corr[j] = crosscorr[j + maxlag];
        ddd++;
    }

    int yn_corr_len = ddd;

    rxx = pmxcorr(noise, noise, p, frameSize);

    for (j = 0; j < p; j++)
    {
        R[j] = rxx[j + p + 1];
        //printf("R[%d], %.32lf\n", j, data->R[j]);

        R_reverse[p - j - 1] = R[j];
        //printf("R_reverse[%d], %.32lf\n", p - j - 1, data->R_reverse[p - j - 1]);
    }

    E0 = rxx[p];

    for (M = 0; M < p; M++)
    {
        if (M == 0)
        {
            K[M] = -R[M] / rxx[p];
            //printf("K[%d], %.32lf\n", M, K[M]);
            a[M] = K[M];
            b[M] = yn_corr[M] / rxx[p];
            //printf("yn_corr[%d], %.32lf\n", M, yn_corr[M]);
            //printf("b[%d], %.32lf\n", M, b[M]);
            //printf("Erxxp = %.32lf\n", rxx[p]);
            E[M] = (1 - (K[M] * K[M]))* rxx[p];
            //printf("K[%d], %.32lf\n", M, K[M]);
            //printf("E[%d], %.32lf\n", M, E[M]);
        }
        else
        {
            temp1 = 0;
            //printf("temp1 = %.32lf\n", temp1);

            for (i = 0; i < M; i++)
            {
                //printf("temp1 = %.32lf\n", temp1);
                temp1 = temp1 + (pre_a[i])* (R_reverse[p - M + i]);
                //printf("temp1 = %.32lf\n", temp1);
                //printf("trial1 = %.32lf\n", pre_a[i] * R_reverse[p - M + i]);

                //printf("pre_a[%d], %.32lf\n", i,  pre_a[i]);
                //printf("R_reverse[%d], %.32lf\n", p - M + i, R_reverse[p - M + i]);
                //printf("temp1 = %.32lf\n", temp1);
                //printf("trial1 = %.32lf\n", pre_a[i] * R_reverse[p - M + i]);
            }

            K[M] = -(temp1 + R[M]) / E[M - 1];
            //printf("R[%d], %.32lf\n", M, R[ M]);
            //printf("K[%d], %.32lf\n", M, K[M]);
            //printf("E[%d], %.32lf\n", M, E[M - 1]);


            a[M] = K[M];

            temp1 = 0;

            temp2 = 0;

            for (i = 0; i < M; i++)
            {
                temp2 = temp2 + pre_b[i] * R_reverse[p - M + i];

                //printf("pre_b[%d], %.32lf\n", i, pre_b[i]);
                //printf(" R_reverse[%d], %.32lf\n", p - M + i, R_reverse[p - M + i]);
                //printf("temp2 = %.32lf\n", temp2);
            }

            b[M] = (yn_corr[M] - temp2) / E[M - 1];
            //printf("b[%d], %.32lf\n", M, b[M]);
            //printf("yn_corr[%d], %.32lf\n", M, yn_corr[M]);
            //printf("E[%d], %.32lf\n", M - 1, E[M - 1]);
            //printf("temp2 = %.32lf\n", temp2);
            temp2 = 0;

            for (i = 0; i < M; i++)
            {
                a[i] = pre_a[i] + K[M] * pre_a[M - 1 - i];

                //printf("pre_a[%d], %.32lf\n", i, pre_a[i]);
                //printf("a[%d], %.32lf\n", i, a[i]);

                b[i] = pre_b[i] + b[M] * pre_a[M - 1 - i];

                //printf("pre_b[%d], %.32lf\n", i, pre_b[i]);
                //printf("b[%d], %.32lf\n", i, b[i]);
            }

            E[M] = (1 - (K[M] * K[M]))*E[M - 1];
            //printf("E[%d], %.32lf\n", M, E[M]);

        }// end of else

        for (i = 0; i < p; i++)
        {
            pre_a[i] = a[i];
            pre_b[i] = b[i];
            //printf("b[%d], %.32lf\n", i, b[i]);
          //  printf("%.32lf\n", b[i]);
            //__android_log_print(ANDROID_LOG_INFO,"bcoff", "bcoff= %1.9g ms",(b[i]));
        }

        //imp[M] =  energy(b,  M);//to find out impulse energy of AFC filter


    }//end of for (M = 0; M < p; M++) loop

    if (E!= NULL)
    {
        free(E);
        E = NULL;
    }
    if (K!= NULL)
    {
        free(K);
        K = NULL;
    }
    if (a!= NULL)
    {
        free(a);
        a = NULL;
    }
    if (pre_a!= NULL)
    {
        free(pre_a);
        pre_a = NULL;
    }
    if (pre_b!= NULL)
    {
        free(pre_b);
        pre_b = NULL;
    }
    if (rxx!= NULL)
    {
        free(rxx);
        rxx = NULL;
    }
    if (R!= NULL)
    {
        free(R);
        R = NULL;
    }
    if (R_reverse!= NULL)
    {
        free(R_reverse);
        R_reverse = NULL;
    }
    if (yn_corr!= NULL)
    {
        free(yn_corr);
        yn_corr = NULL;
    }
    if (imp!= NULL)
    {
        free(imp);
        imp= NULL;
    }


    return b;

}//end of GLD function

double multi(double* x, double *y, int length)
{
    double sum = 0;

    for (int i = 0; i < length; i++)
    {
        sum = sum + x[i] * y[i];

    }
    return sum;
}

//function for implementing feedback canceller
double canceller(double mic_in, double *bcoff, int maxlag, double *V_f)
{
    int filter_length = maxlag + 1;

    double z = 0;
    float a=0;
    float b=0;
    //float pp[200];


    for (int j = filter_length - 1; j >= 0; j--)
    {
        if (j == 0)
        {
            V_f[j] = mic_in;
           // pp[j] = V_f[j];
        }
        else
        {
            V_f[j] = V_f[j - 1];
           // pp[j] = V_f[j];
        }
    }

    z = multi(bcoff, V_f, filter_length);
/*
    for (int jj =0; jj < maxlag + 1+1; jj++)
    {
        if (jj == 0)
        {
            Buffer[jj] = z;
        }
        else
        {
            Buffer[jj] = V_f[jj];
        }

    }
*/
    return z;
}

double *canceller2(double mic_in, int maxlag, double *V_f)
{
    int filter_length = maxlag + 1;
    float a=0;
    float b=0;
    float pp2[200];

    for (int j = filter_length - 1; j >= 0; j--)
    {
        if (j == 0)
        {
            V_f[j] = mic_in;
            pp2[j] = V_f[j];
        }
        else
        {
            V_f[j] = V_f[j - 1];
            pp2[j] = V_f[j];
        }
    }

    //z = multi(bcoff, V_f, filter_length);

    return V_f;
}