// FIR_Hamming_Filter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <vector> 
#include <complex>
#include <algorithm>
#include "FIR_Hamming_Filter.h"

#define PI 3.141592653589793

int main()
{
    
    double f1 = 1;  //High Pass
    double f2 = 20; //Low Pass
    double sf = 2048;    //Sampling frequency
    int order_filt = 200; //Order
    
    std::vector<std::vector<double> > coeff_final(2);
   
    int type_filt = 3;
    FIR_H_F::FIR_Hamming_Filter fir;
    
    switch (type_filt)
    {
    case 0:


        coeff_final = fir.Fir_LP(order_filt, sf, f1);
       
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

            for (int ll = 0; ll < coeff_final[0].size(); ll++)

            {
                std::cout << coeff_final[kk][ll] << std::ends;

            }

            std::cout << std::endl;

        }

        break;

    case 1:

        coeff_final = fir.Fir_HP(order_filt, sf, f1);

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

            for (int ll = 0; ll < coeff_final[0].size(); ll++)

            {

                std::cout << coeff_final[kk][ll] << std::ends;

            }

            std::cout << std::endl;

        }

        break;

    case 2:

        coeff_final = fir.Fir_BP(order_filt, sf, f1, f2);
       
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

            for (int ll = 0; ll < coeff_final[0].size(); ll++)

            {

                std::cout << coeff_final[kk][ll] << std::ends;

            }

            std::cout << std::endl;
        }

        break;

    case 3:
       
        coeff_final = fir.Fir_BS(order_filt, sf, f1, f2);

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

            for (int ll = 0; ll < coeff_final[0].size(); ll++)

            {

                std::cout << coeff_final[kk][ll] << std::ends;

            }

            std::cout << std::endl;

        }

        break;


    }

    int length_signal = 2048;
    std::vector<double> signal;
    for (int kk = 0; kk < sf; kk++)
    {

        signal.push_back(0);

    }
    

    for (int kk = 0; kk < length_signal; kk++)
    {

        signal[kk] = std::sin(2 * PI * 10 * (double)kk / (double)sf) + std::sin(2 * PI * 45 * (double)kk / (double)sf);

    }

    std::vector<double> filt_sign = fir.Filter_Data(coeff_final, signal);

    for (int kk = 0; kk < length_signal; kk++)
    {

        std::cout << filt_sign[kk] << "\n";

    }

}

