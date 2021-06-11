#pragma once
#include <stdio.h> 
#include <complex>
#include <vector>

#ifndef FIR_Hamming_Filter_H
#define FIR_Hamming_Filter_H

#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif

    namespace FIR_H_F
    {

        class FIR_Hamming_Filter

        {
            

        private:

            // Compute the vector of frequencies to pass to FIRLS
            void Desired_Freq(double);
            void Desired_Freq(double, double);

            // Compute the magnitude vector
            void Desired_MAG(int, int);

            // Build a Hamming Window
            void Hamm_win(int);

            // Sinc function
            double Sinc_Function(double);

            // Compute windowed impulse response
            void firls(int, int, std::vector<double>, std::vector<double>);

            // Scale the filter coefficients
            void Scale_Filter();

        public:

            //Estimate the coeffients of a band-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator. Order, sf, f1, f2
            std::vector<std::vector<double> > Fir_BP(int, double, double, double);

            //Estimate the coeffients of a band-stop filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator. Order, sf, f1, f2
            std::vector<std::vector<double> > Fir_BS(int, double, double, double);

            //Estimate the coeffients of a low-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator. Order, sf, f1
            std::vector<std::vector<double> > Fir_LP(int, double, double);

            //Estimate the coeffients of a high-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator. Order, sf, f1
            std::vector<std::vector<double> > Fir_HP(int, double, double);

            //Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
            std::vector<double> Filter_Data(std::vector<std::vector<double> > coeff_filt, std::vector<double> pre_filt_signal);

        };

    }

#endif

#ifdef __cplusplus

}

#endif
