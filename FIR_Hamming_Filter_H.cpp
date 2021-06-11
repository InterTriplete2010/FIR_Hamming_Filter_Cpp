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

//Global variables
int nbands;
int first_band;
std::vector<double> freq;
std::vector<double> mags;
std::vector<double> magnitude;
std::vector<double> hamm_window;
std::vector<double> F_Matlab;
std::vector<double> h_Matlab;
std::vector<std::vector<double>> b_filter;


using namespace FIR_H_F;

//-------------------------------------------------------------------------------------------------------------------------------------//
 //Step 1: Compute the vector of frequencies to pass to FIRLS
void FIR_Hamming_Filter::Desired_Freq(double Wnf_1)
{

    nbands = 2;

    //Initialize the vector freq
    for (int i = 0; i < nbands * 2; i++)
    {

        freq.push_back(0);

    }

    //create the pair of frequency points based on Wn
    freq[0] = 0;
    freq[1] = Wnf_1;
    freq[2] = Wnf_1;
    freq[3] = 1;
}

//Overload of Step 1
void FIR_Hamming_Filter::Desired_Freq(double Wnf_1, double Wnf_2)
{

    nbands = 3;

    //Initialize the vector freq
    for (int i = 0; i < nbands * 2; i++)
    {

        freq.push_back(0);
   
    }

    //create the pair of frequency points based on Wn
    freq[0] = 0;
    freq[1] = Wnf_1;
    freq[2] = Wnf_1;
    freq[3] = Wnf_2;
    freq[4] = Wnf_2;
    freq[5] = 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------------//
//Step 2: Compute the magnitude vector
void FIR_Hamming_Filter::Desired_MAG(int type_filt, int nbands)
{

    if (type_filt == 1 || type_filt == 2)
    {

        first_band = 0;

    }

    else
    {

        first_band = 1;

    }

    //Initialize the vector mags
    for (int i = 0; i < nbands; i++)
    {

        mags.push_back(0);

    }
    

    for (int kk = 0; kk < nbands; kk++)
    {

        mags[kk] = (first_band + kk) % 2;

    }

    int track_mags = 0;
    //Initialize the vector magnitude
    for (int i = 0; i < nbands * 2; i++)
    {

        magnitude.push_back(0);

    }
    

    for (int kk = 0; kk < nbands * 2; kk += 2)
    {

        magnitude[kk] = mags[track_mags];
        magnitude[kk + 1] = mags[track_mags];

        track_mags++;

    }

}

//Step 3: Build a Hamming Window
void FIR_Hamming_Filter::Hamm_win(int filter_order)
{

    double coeff_I = 0.54;
    double coeff_II = 0.46;

    //Initialize the vector hamm_window
    for (int i = 0; i < filter_order; i++)
    {

        hamm_window.push_back(0);

    }
    

    for (int kk = 0; kk < filter_order; kk++)
    {

        hamm_window[kk] = coeff_I - coeff_II * std::cos(2 * PI * kk / (filter_order - 1));

    }

}
//-------------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------------//
//Step 4: Sinc function
double FIR_Hamming_Filter::Sinc_Function(double input_value)
{
    double output_sinc;

    if (input_value == 0)
    {

        output_sinc = 1;

    }

    else
    {

        output_sinc = std::sin(input_value * PI) / (input_value * PI);

    }

    return output_sinc;

}
//-------------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------------//
//Step 5: Compute windowed impulse response
void FIR_Hamming_Filter::firls(int filter_order, int filter_type, std::vector<double> freq, std::vector<double> magnitude)
{
    bool odd_order = false;

    if ((filter_order % 2) != 0)
    {

        filter_order += 1;
        odd_order = true;

    }

    double freq_length = freq.size();
    double length_weight = std::floor(freq_length / 2);
    std::vector<double> weight;
    
    //Initialize the vector weight
    for( int kk = 0; kk < length_weight; kk++)
    {

        weight.push_back(0);

    }


    for (int kk = 0; kk < weight.size(); kk++)
    {

        weight[kk] = 1;

    }

    //Initialize the vector F_Matlab
    for (int kk = 0; kk < freq.size(); kk++)
    {

        F_Matlab.push_back(0);

    }
    
    for (int kk = 0; kk < F_Matlab.size(); kk++)
    {

        F_Matlab[kk] = freq[kk] / 2;

    }

    std::vector<double> amp_Matlab;
    
    //Initialize the vector amp_Matlab
    for (int kk = 0; kk < magnitude.size(); kk++)
    {

        amp_Matlab.push_back(0);

    }
   
    amp_Matlab = magnitude;

    std::vector<double> wt_Matlab;
    
    //Initialize the vector wt_Matlab
    for (int kk = 0; kk < length_weight; kk++)
    {

        wt_Matlab.push_back(0);

    }

   
    for (int kk = 0; kk < wt_Matlab.size(); kk++)
    {

        wt_Matlab[kk] = std::abs(std::sqrt(weight[kk]));

    }

    //Calculate the "diff" of F_Matlab
    std::vector<double> dF_Matlab;

    //Initialize the vector dF_Matlab
    for (int kk = 0; kk < freq.size() - 1; kk++)
    {

        dF_Matlab.push_back(0);

    }

    for (int kk = 0; kk < dF_Matlab.size(); kk++)
    {

        dF_Matlab[kk] = F_Matlab[kk + 1] - F_Matlab[kk];

    }

    int lendF = freq.size() - 1;

    // validate weight
    std::vector<double> temp_W;

    //Initialize the vector temp_W
    for (int kk = 0; kk < wt_Matlab.size() - 1; kk++)
    {

        temp_W.push_back(0);

    }

    for (int kk = 0; kk < temp_W.size(); kk++)
    {

        temp_W[kk] = wt_Matlab[kk] - wt_Matlab[0];

    }

    // find the order
    int L_Matlab = filter_order / 2;

    // odd order
    bool Nodd = false;

    if (odd_order == false)
    {

        Nodd = (((filter_order + 1) % 2) == 1);

    }

    else
    {

        Nodd = (((filter_order) % 2) == 1);

    }

    // initialize b0
    double b0_Matlab = 0;

    filter_type = 0;

    // Type I and Type II linear phase FIR
    std::vector<double> m_Matlab;

    if (odd_order == true)
    {
            //Initialize the vector m_Matlab
            for (int kk = 0; kk < L_Matlab; kk++)
            {

                m_Matlab.push_back(0);

            }

    }

    else
    {

        //Initialize the vector m_Matlab
        for (int kk = 0; kk < L_Matlab + 1; kk++)
        {

            m_Matlab.push_back(0);

        }

    }

    if (filter_type == 0)
    {
        // Basis vectors are cos(2 * pi * m * f)
        if (!Nodd == true)
        {

            for (int kk = 0; kk < m_Matlab.size(); kk++)
            {

                m_Matlab[kk] = kk + 0.5;

            }

        }

        else
        {

            for (int kk = 0; kk < m_Matlab.size(); kk++)
            {

                m_Matlab[kk] = kk;

            }

        }

    }

    std::vector<double> k_Matlab;
   
    if (odd_order == true)
    {

        for (int kk = 0; kk < m_Matlab.size(); kk++)
        {

            k_Matlab.push_back(0);

        }
      
    }

    else
    {

        for (int kk = 0; kk < m_Matlab.size() - 1; kk++)
        {

            k_Matlab.push_back(0);

        }

    }


    if (Nodd == true)
    {

        for (int kk = 0; kk < k_Matlab.size(); kk++)
        {

            k_Matlab[kk] = m_Matlab[kk + 1];

        }

        b0_Matlab = 0; //First entry must be handled separately(where k(1) = 0)

    }

    else
    {

        k_Matlab = m_Matlab;

    }

    // preallocate b matrix
    std::vector<double> b_Matlab;
    for (int kk = 0; kk < k_Matlab.size(); kk++)
    {

        b_Matlab.push_back(0);

    }


    double m_s_Matlab;
    double b1_Matlab;
    
    std::vector<double> A_Matlab;
    for (int kk = 0; kk < magnitude.size(); kk++)
    {

        A_Matlab.push_back(0);

    }
    
    A_Matlab = amp_Matlab;

    for (int kk = 0; kk < F_Matlab.size(); kk += 2)
    {

        m_s_Matlab = (A_Matlab[kk + 1] - A_Matlab[kk]) / (F_Matlab[kk + 1] - F_Matlab[kk]);   // Slope
        b1_Matlab = A_Matlab[kk + 1] - m_s_Matlab * F_Matlab[kk];                           // y - intercept

        if (Nodd == true)
        {

            b0_Matlab = b0_Matlab + (b1_Matlab * (F_Matlab[kk + 1] - F_Matlab[kk]) + m_s_Matlab / 2 * (F_Matlab[kk + 1] * F_Matlab[kk + 1] - F_Matlab[kk] * F_Matlab[kk])) * (std::pow((wt_Matlab[(kk + 2) / 2 - 1]), 2));

        }

        for (int ll = 0; ll < b_Matlab.size(); ll++)
        {

            b_Matlab[ll] += (m_s_Matlab / (4 * PI * PI) * (std::cos(2 * PI * k_Matlab[ll] * F_Matlab[kk + 1]) - std::cos(2 * PI * k_Matlab[ll] * F_Matlab[kk])) / (k_Matlab[ll] * k_Matlab[ll])) * (std::pow(wt_Matlab[(kk + 2) / 2 - 1], 2));
            b_Matlab[ll] += (F_Matlab[kk + 1] * (m_s_Matlab * F_Matlab[kk + 1] + b1_Matlab) * Sinc_Function(2 * k_Matlab[ll] * F_Matlab[kk + 1]) - F_Matlab[kk] * (m_s_Matlab * F_Matlab[kk] + b1_Matlab) * Sinc_Function(2 * k_Matlab[ll] * F_Matlab[kk])) * (std::pow(wt_Matlab[(kk + 2) / 2 - 1], 2));
        }

    }

    //Increase the size of the array and add b0_Matlab in the first cell if Nodd is true
    if (Nodd == true)
    {
        b_Matlab.push_back(0);

        for (int hh = 0; hh < b_Matlab.size(); hh++)
        {

            if (hh == b_Matlab.size() - 1)
            {

                b_Matlab[b_Matlab.size() - hh - 1] = b0_Matlab;

            }

            else
            {

                b_Matlab[b_Matlab.size() - hh - 1] = b_Matlab[b_Matlab.size() - hh - 2];

            }

        }
    }

    std::vector<double> a_Matlab;
    for (int kk = 0; kk < b_Matlab.size(); kk++)
    {

        a_Matlab.push_back(0);

    }
    
    for (int hh = 0; hh < a_Matlab.size(); hh++)
    {

        a_Matlab[hh] = std::pow(wt_Matlab[0], 2) * 4 * b_Matlab[hh];

        if (hh == 0 && Nodd == true)
        {

            a_Matlab[hh] /= 2;

        }

    }

    //-------------------------------------------------------------------------------------//
    //Compute the unwindowed impulse response of the filter
    if (Nodd == true)
    {

        for (int kk = 0; kk < filter_order + 1; kk++)
        {

            h_Matlab.push_back(0);

        }
        
        for (int hh = 0; hh < h_Matlab.size(); hh++)
        {

            if (hh < filter_order / 2)
            {

                h_Matlab[hh] = a_Matlab[a_Matlab.size() - hh - 1] / 2;


            }

            if (hh == filter_order / 2)
            {

                h_Matlab[hh] = a_Matlab[0];

            }

            if (hh > filter_order / 2)
            {

                h_Matlab[hh] = a_Matlab[hh - h_Matlab.size() / 2] / 2;

            }

        }

        //-------------------------------------------------------------------------------------//
        //Compute the windowed (Hamming) impulse response of the filter
        Hamm_win(h_Matlab.size());
        for (int kk = 0; kk < h_Matlab.size(); kk++)
        {

            h_Matlab[kk] = h_Matlab[kk] * hamm_window[kk];

        }
        //-------------------------------------------------------------------------------------//

    }

    else
    {

        if (odd_order == true)
        {

            for (int kk = 0; kk < filter_order; kk++)
            {

                h_Matlab.push_back(0);

            }
            
        }

        else
        {

            for (int kk = 0; kk < filter_order + 1; kk++)
            {

                h_Matlab.push_back(0);

            }
           
        }

        for (int hh = 0; hh < h_Matlab.size(); hh++)
        {

            if (hh < h_Matlab.size() / 2)
            {

                h_Matlab[hh] = a_Matlab[a_Matlab.size() - hh - 1] * 0.5;

            }

            else
            {

                h_Matlab[hh] = a_Matlab[hh - a_Matlab.size()] * 0.5;

            }

        }

        //-------------------------------------------------------------------------------------//
        //Compute the windowed (Hamming) impulse response of the filter
        Hamm_win(h_Matlab.size());
        for (int kk = 0; kk < h_Matlab.size(); kk++)
        {

            h_Matlab[kk] = h_Matlab[kk] * hamm_window[kk];

        }
        //-------------------------------------------------------------------------------------//

    }

}

//-------------------------------------------------------------------------------------//
//Step 6: Scale the filter coefficients
void FIR_Hamming_Filter::Scale_Filter()
{
    double f0_Matlab;
    //b_filter = new double[h_Matlab.Length];

    //Initialize the matrix where to save the coefficients
    std::vector<double> temp_v;

    for (int ff = 0; ff < h_Matlab.size(); ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        b_filter.push_back(temp_v);

    }
    

    //Setting the first value of the numerator to "1"
    b_filter[1][0] = 1;

    if (first_band == 1)
    {

        double temp_sum = 0;

        for (int hh = 0; hh < h_Matlab.size(); hh++)
        {

            temp_sum += h_Matlab[hh];

        }

        for (int hh = 0; hh < h_Matlab.size(); hh++)
        {

            b_filter[0][hh] = h_Matlab[hh] / temp_sum;

        }

    }

    else
    {

        if (freq[3] == 1)
        {

            // Unity gain at Fs / 2
            f0_Matlab = 1;

        }

        else
        {

            // unity gain at center of first passband
            f0_Matlab = (freq[2] + freq[3]) / 2;

        }
        
        std::complex<double> complex_imag (0.0, 1.0);
        std::complex<double> temp_den (0.0, 0.0);
        for (int kk = 0; kk < h_Matlab.size(); kk++)
        {

            temp_den += std::exp(-complex_imag * 2.0 * PI * (double)kk * (f0_Matlab / 2)) * h_Matlab[kk];

        }

        double temp_den_double = std::abs(temp_den);

        for (int kk = 0; kk < h_Matlab.size(); kk++)
        {

            b_filter[0][kk] = h_Matlab[kk] / temp_den_double;

        }
    }
}
//-------------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------------//
//Calculate the coefficients
//Low-Pass filter
std::vector<std::vector<double>> FIR_Hamming_Filter::Fir_LP(int filt_order, double sf, double f1)
{

    //Clean up the global variables for a new analysis
    if (b_filter.size() > 0)
    {

        freq.erase(freq.begin(), freq.begin() + freq.size());
        mags.erase(mags.begin(), mags.begin() + freq.size());
        magnitude.erase(magnitude.begin(), magnitude.begin() + magnitude.size());
        hamm_window.erase(hamm_window.begin(), hamm_window.begin() + hamm_window.size());
        F_Matlab.erase(F_Matlab.begin(), F_Matlab.begin() + F_Matlab.size());
        h_Matlab.erase(h_Matlab.begin(), h_Matlab.begin() + h_Matlab.size());
        
        b_filter.erase(b_filter.begin(), b_filter.begin() + b_filter.size());
       
    }

    int filter_type = 0;
    //Normalizing the corner frequency with respect to the nyquist frequency
    f1 = f1 / (sf / 2);

    //Check that the normalized frequencies are within the correct range of values
    if (f1 <= 0 || f1 >= 1)
    {

        throw ("Cut-off frequencies must be in the (0,1) range");

    }

    //Check that the order of the filter is > 0
    if (filt_order <= 0)
    {

        throw ("The order of the filter must be > 0");

    }

    Desired_Freq(f1);
    Desired_MAG(filter_type, nbands);
    firls(filt_order, filter_type, freq, magnitude);
    Scale_Filter();

    return b_filter;

}

//High-Pass filter
std::vector<std::vector<double>> FIR_Hamming_Filter::Fir_HP(int filt_order, double sf, double f1)
{

    //Clean up the global variables for a new analysis
    if (b_filter.size() > 0)
    {

        freq.erase(freq.begin(), freq.begin() + freq.size());
        mags.erase(mags.begin(), mags.begin() + freq.size());
        magnitude.erase(magnitude.begin(), magnitude.begin() + magnitude.size());
        hamm_window.erase(hamm_window.begin(), hamm_window.begin() + hamm_window.size());
        F_Matlab.erase(F_Matlab.begin(), F_Matlab.begin() + F_Matlab.size());
        h_Matlab.erase(h_Matlab.begin(), h_Matlab.begin() + h_Matlab.size());

        b_filter.erase(b_filter.begin(), b_filter.begin() + b_filter.size());

    }

    int filter_type = 1;

    int temp_odd = filt_order % 2;
    if (temp_odd == 1)
    {

        //Increasing the order of the filter by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency
        filt_order += 1;

    }

    //Normalizing the corner frequency with respect to the nyquist frequency
    f1 = f1 / (sf / 2);

    //Check that the normalized frequencies are within the correct range of values
    if (f1 <= 0 || f1 >= 1)
    {

        throw ("Cut-off frequencies must be in the (0,1) range");

    }

    //Check that the order of the filter is > 0
    if (filt_order <= 0)
    {

        throw ("The order of the filter must be > 0");

    }

    Desired_Freq(f1);
    Desired_MAG(filter_type, nbands);
    firls(filt_order, filter_type, freq, magnitude);
    Scale_Filter();

    return b_filter;

}

//Band-Pass filter
std::vector<std::vector<double>> FIR_Hamming_Filter::Fir_BP(int filt_order, double sf, double f1, double f2)
{

    //Clean up the global variables for a new analysis
    if (b_filter.size() > 0)
    {

        freq.erase(freq.begin(), freq.begin() + freq.size());
        mags.erase(mags.begin(), mags.begin() + freq.size());
        magnitude.erase(magnitude.begin(), magnitude.begin() + magnitude.size());
        hamm_window.erase(hamm_window.begin(), hamm_window.begin() + hamm_window.size());
        F_Matlab.erase(F_Matlab.begin(), F_Matlab.begin() + F_Matlab.size());
        h_Matlab.erase(h_Matlab.begin(), h_Matlab.begin() + h_Matlab.size());

        b_filter.erase(b_filter.begin(), b_filter.begin() + b_filter.size());

    }

    int filter_type = 2;

    //Normalizing the corner frequencies with respect to the nyquist frequency
    f1 = f1 / (sf / 2);
    f2 = f2 / (sf / 2);

    //Check that the normalized frequencies are within the correct range of values
    if (f1 <= 0 || f1 >= 1 || f2 <= 0 || f2 >= 1)
    {

        throw ("Cut-off frequencies must be in the (0,1) range");

    }

    //Check that f1 < f2
    if (f1 > f2)
    {

        throw ("The corner frequency f1 needs to be lower than f2");

    }

    //Check that the order of the filter is > 0
    if (filt_order <= 0)
    {

        throw ("The order of the filter must be > 0");

    }

    Desired_Freq(f1, f2);
    Desired_MAG(filter_type, nbands);
    firls(filt_order, filter_type, freq, magnitude);
    Scale_Filter();

    return b_filter;

}

//Band-Stop filter
std::vector<std::vector<double>> FIR_Hamming_Filter::Fir_BS(int filt_order, double sf, double f1, double f2)
{

    //Clean up the global variables for a new analysis
    if (b_filter.size() > 0)
    {

        freq.erase(freq.begin(), freq.begin() + freq.size());
        mags.erase(mags.begin(), mags.begin() + freq.size());
        magnitude.erase(magnitude.begin(), magnitude.begin() + magnitude.size());
        hamm_window.erase(hamm_window.begin(), hamm_window.begin() + hamm_window.size());
        F_Matlab.erase(F_Matlab.begin(), F_Matlab.begin() + F_Matlab.size());
        h_Matlab.erase(h_Matlab.begin(), h_Matlab.begin() + h_Matlab.size());

        b_filter.erase(b_filter.begin(), b_filter.begin() + b_filter.size());

    }

    int filter_type = 3;

    int temp_odd = filt_order % 2;
    if (temp_odd == 1)
    {

        //Increasing the order of the filter by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency
        filt_order += 1;

    }

    //Normalizing the corner frequencies with respect to the nyquist frequency
    f1 = f1 / (sf / 2);
    f2 = f2 / (sf / 2);

    //Check that the normalized frequencies are within the correct range of values
    if (f1 <= 0 || f1 >= 1 || f2 <= 0 || f2 >= 1)
    {

        throw ("Cut-off frequencies must be in the (0,1) range");

    }

    //Check that f1 < f2
    if (f1 > f2)
    {

        throw ("The corner frequency f1 needs to be lower than f2");

    }

    //Check that the order of the filter is > 0
    if (filt_order <= 0)
    {

        throw ("The order of the filter must be > 0");

    }

    Desired_Freq(f1, f2);
    Desired_MAG(filter_type, nbands);
    firls(filt_order, filter_type, freq, magnitude);
    Scale_Filter();

    return b_filter;

}


//-------------------------------------------------------------------------------------------------------------------------------------//
//Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
std::vector<double> FIR_Hamming_Filter::Filter_Data(std::vector<std::vector<double> > coeff_filt, std::vector<double> pre_filt_signal)
{

    std::vector<double> filt_signal(pre_filt_signal.size(), 0.0);

    std::vector<std::vector<double>> w_val;
    std::vector<double> temp_v;

    for (int ff = 0; ff < pre_filt_signal.size(); ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < coeff_filt[0].size(); hh++)
    {

        w_val.push_back(temp_v);

    }


    //Convolution product to filter the data
    for (int kk = 0; kk < pre_filt_signal.size(); kk++)
    {

        if (kk == 0)
        {

            filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0];

            for (int ww = 1; ww < coeff_filt[0].size(); ww++)
            {

                w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];

            }

        }

        else
        {

            filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0] + w_val[0][kk - 1];

            for (int ww = 1; ww < coeff_filt[0].size(); ww++)
            {

                w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] + w_val[ww][kk - 1] - filt_signal[kk] * coeff_filt[1][ww];

                if (ww == coeff_filt[0].size() - 1)
                {

                    w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];

                }

            }

        }

    }


    return filt_signal;

}

//-------------------------------------------------------------------------------------------------------------------------------------//