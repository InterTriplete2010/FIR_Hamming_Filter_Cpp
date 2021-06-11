# FIR_Hamming_Filter_Cpp
FIR filter based on a Hamming window

C# Library to calculate the coefficients of a FIR filter with Hamming window and to filter the data

This code calculates the coefficients of the Band-pass, Band-stop, Low-pass and High-pass FIR filters. It also filters the data, but no zero-phase delay is applied. The namespace is: FIR_H_F.

Each filter function will return a 2 rows x N coefficients 2D vector, where Row 1 = Numerator and Row 2 = Denumerator.

Band-pass: the function is "std::vector<std::vector<double> > Fir_BP(int, double, double, double)". The first argument if the order of the filter, the second one if the sampling frequency, while the last two arguments are the two cut-off frequencies f1 and f2;

Band-stop: the function is "std::vector<std::vector<double> > Fir_BS(int, double, double, double)". The first argument if the order of the filter, the second one if the sampling frequency, while the last two arguments are the two cut-off frequencies f1 and f2. If the order of the filter is odd, it will be increased by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency;

Low-pass: the function is "std::vector<std::vector<double> > Fir_LP(int, double, double)". The first argument if the order of the filter, the second one if the sampling frequency, while the last one is the cut-off frequency f1;

High-pass: the function is "std::vector<std::vector<double> > Fir_HP(int, double, double)". The first argument if the order of the filter, the second one if the sampling frequency, while the last one is the cut-off frequency f1. If the order of the filter is odd, it will be increased by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency;

Filter the data: the method is "std::vector<double> Filter_Data(std::vector<std::vector<double> > coeff_filt, std::vector<double> pre_filt_signal)". The two arguments are the filter coefficients and the signal to be filtered. It returns the filtered signal.

If you have any question and/or want to report bugs, please e-mail me (Ale) at: pressalex@hotmail.com
