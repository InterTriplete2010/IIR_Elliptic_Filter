//This file contains the 'main' function. Program execution begins and ends there.
//Calculate the coefficients of the butterworth filter. It follows the steps as described in Matlab. 
//The same symbols as in Matlab are used to facilitate the writing of the code
#define _CRTDBG_MAP_ALLOC
#include <stdio.h> 
#include <stdlib.h> 
#include <crtdbg.h>
#include <iostream>
#include <math.h>
#include <string.h> 
#include <vector> 
#include <complex.h>
#include "Elliptic_Filter.h"

#define ARMA_DONT_USE_CXX11
#include <armadillo>

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

#define PI 3.141592653589793

int main()
{


    double f1 = 20;  //High Pass
    double f2 = 40; //Low Pass
    double Rp = 8;  //Decibels of peak - to - peak passband ripple 
    double Rs = 30; //Decibels of stopband attenuation down from the peak passband value
    double sf = 2048;    //Sampling frequency
    int order_filt = 3; //Order
    double Nyquist_F = sf / 2;

    std::vector<std::vector<double> > coeff_final(2);

    int type_filt = 0;
    IIR_E_F::IIR_Elliptic_Filter ir_b;

    bool check_stability_flag;

    //_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
    //_CrtSetBreakAlloc(154);
    //_CrtSetBreakAlloc(155);
    //_CrtSetBreakAlloc(156);

    int coeff_numb = 0;
    switch (type_filt)
    {
    case 0:

        coeff_numb = 2 * (size_t)order_filt + 1;


        for (int i = 0; i < 2; i++)
        {

            coeff_final[i].resize(coeff_numb);

        }


        _CrtDumpMemoryLeaks();

        coeff_final = ir_b.lp2bp(order_filt, Rp, Rs, f1 / Nyquist_F, f2 / Nyquist_F);

        check_stability_flag = ir_b.check_stability_iir(coeff_final);

        if (check_stability_flag)
        {

            std::cout << "The filter is stable" << std::endl;

        }

        else
        {

            std::cout << "The filter is unstable" << std::endl;

        }

        for (int kk = 0; kk < 2; kk++)
        {
            if (kk == 0)
            {

                std::cout << "Numerator: " << std::ends;

            }

            else
            {

                std::cout << "Denumerator: " << std::ends;

            }

            for (int ll = 0; ll < coeff_numb; ll++)

            {
                std::cout << coeff_final[kk][ll] << "\t" << std::ends;

            }

            std::cout << std::endl;

        }

        break;

    case 1:

        coeff_numb = 2 * order_filt + 1;


        for (int i = 0; i < 2; i++)
        {

            coeff_final[i].resize(coeff_numb);

        }

        _CrtDumpMemoryLeaks();
        coeff_final = ir_b.lp2bs(order_filt, Rp, Rs, f1 / Nyquist_F, f2 / Nyquist_F);


        check_stability_flag = ir_b.check_stability_iir(coeff_final);

        if (check_stability_flag)
        {

            std::cout << "The filter is stable" << std::endl;

        }

        else
        {

            std::cout << "The filter is unstable" << std::endl;

        }

        for (int kk = 0; kk < 2; kk++)
        {

            if (kk == 0)
            {

                std::cout << "Numerator: " << std::ends;

            }

            else
            {

                std::cout << "Denumerator: " << std::ends;

            }

            for (int ll = 0; ll < coeff_numb; ll++)

            {

                std::cout << coeff_final[kk][ll] << "\t" << std::ends;

            }

            std::cout << std::endl;

        }

        break;

    case 2:

        coeff_numb = order_filt + 1;

        for (int i = 0; i < 2; i++)
        {

            coeff_final[i].resize(coeff_numb);

        }

        _CrtDumpMemoryLeaks();
        coeff_final = ir_b.lp2hp(order_filt, Rp, Rs, f1 / Nyquist_F);

        check_stability_flag = ir_b.check_stability_iir(coeff_final);

        if (check_stability_flag)
        {

            std::cout << "The filter is stable" << std::endl;

        }

        else
        {

            std::cout << "The filter is unstable" << std::endl;

        }

        for (int kk = 0; kk < 2; kk++)
        {

            if (kk == 0)
            {

                std::cout << "Numerator: " << std::ends;

            }

            else
            {

                std::cout << "Denumerator: " << std::ends;

            }

            for (int ll = 0; ll < coeff_numb; ll++)

            {

                std::cout << coeff_final[kk][ll] << "\t" << std::ends;

            }

            std::cout << std::endl;
        }

        break;

    case 3:
        coeff_numb = order_filt + 1;

        for (int i = 0; i < 2; i++)
        {

            coeff_final[i].resize(coeff_numb);

        }

        _CrtDumpMemoryLeaks();
        coeff_final = ir_b.lp2lp(order_filt, Rp, Rs, f2 / Nyquist_F);

        check_stability_flag = ir_b.check_stability_iir(coeff_final);

        if (check_stability_flag)
        {

            std::cout << "The filter is stable" << std::endl;

        }

        else
        {

            std::cout << "The filter is unstable" << std::endl;

        }

        for (int kk = 0; kk < 2; kk++)
        {

            if (kk == 0)
            {

                std::cout << "Numerator: " << std::ends;

            }

            else
            {

                std::cout << "Denumerator: " << std::ends;

            }

            for (int ll = 0; ll < coeff_numb; ll++)

            {

                std::cout << coeff_final[kk][ll] << "\t" << std::ends;

            }

            std::cout << std::endl;

        }

        break;


    }

    //Create a complex sine wave
    std::vector<double> test_sin(sf, 0.0);
    double sf_1 = 20;
    double sf_2 = 60;

    //Write the original waveform to a file
    std::ofstream myfile;
    myfile.open("Waveform.txt");
    
    for (double kk = 0; kk < sf; kk++)
    {

        test_sin[kk] = sin(2 * PI * kk * sf_1 / sf) + sin(2 * PI * kk * sf_2 / sf);

        myfile << test_sin[kk];
        myfile << "\n";

    }

    std::vector<double> filt_sign = ir_b.Filter_Data(coeff_final, test_sin);

    //Write the output of the filter to a file
    std::ofstream myfile_I;
    myfile_I.open("Output_Filter.txt");
    for (int kk = 0; kk < filt_sign.size(); kk++)
    {

        myfile_I << filt_sign[kk];
        myfile_I << "\n";

    }

}
