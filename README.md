# IIR_Elliptic_Filter
C++ code to calculate the coefficients of the Elliptic filter


This code calculates the coefficients of the Band-pass, Band-stop, Low-pass and High-pass Elliptic filters. The file Elliptic_Filter.cpp can be used to test the code. It also filters the data, but no zero-phase delay is applied. The name space is: IIR_E_F. The code follows the same steps as in Matlab.

Each filter function will return a 2 rows x N coefficients 2D vector, where Row 1 = Numerator and Row 2 = Denumerator. The method "check_stability_iir" can be used to check the stability of the filter.

Band-pass: the function is "std::vector<std::vector > lp2bs(int, double, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last two arguments are the two normalized cut-off frequencies (f1/NF, f2/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1). Please, keep in mind that if you enter order_filt = 2, the order of the filter will be 2 * order_filt = 4;

Band-stop: the function is "std::vector<std::vector > lp2bs(int, double, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last two arguments are the two normalized cut-off frequencies (f1/NF, f2/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1). Please, keep in mind that if you enter order_filt = 2, the order of the filter will be 2 * order_filt = 4;

High-pass: the function is "std::vector<std::vector > lp2hp(int, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last argument is the normalized cut-off frequency (f/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1);

Low-pass: the function is "std::vector<std::vector > lp2hp(int, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last argument is the normalized cut-off frequency (f/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1);

Check the stability of the filter: the method is "bool check_stability_iir(std::vector<std::vector >)". The argument is the 2D array containing the filter coefficients. It returns "true" if the filter is stable, "false" if it is unstable.

Filter the data: the method is "std::vector Filter_Data(std::vector<std::vector > coeff_filt, std::vector pre_filt_signal)". The two arguments are the filter coefficients and the signal to be filtered. It returns the filtered signal.

The library Armadillo needs to be downloaded and installed ((http://arma.sourceforge.net/download.html)) along with lapack and blas. I have uploaded the lapack and blas libraries that I have used. Please, note that with the older version of Armadillo, I had to use #define ARMA_DONT_USE_CXX11 to make armadillo library work with C++/CLR in visual studio 2019. If you use the latest version (armadillo-9.880.1), which I would recommend, because it is supposedly faster than the previous one, as the developers told me, you should replace #define ARMA_DONT_USE_CXX11 with #define ARMA_DONT_USE_CXX11_MUTEX.

If you are running the code in Linux (I did it in Ubuntu 18.04), you need to make the following changes in Elliptic_Filter.cpp:

Comment out #include <crtdbg.h>
Cooment out _CrtDumpMemoryLeaks() in lines 68, 129, 191 and 250.
You can also remove the line of code #define ARMA_DONT_USE_CXX11_MUTEX (or #define ARMA_DONT_USE_CXX11) and, based on my communication with the guys who created armadillo, the code should run faster.

Compile the code in the following way:

g++ -c Elliptic_Filter_H.cpp (to generate Elliptic_Filter_H.o)
g++ -ggdb Elliptic_Filter.cpp Elliptic_Filter_H.o -larmadillo -o <Name_Exe_File>

If you have any questions and/or want to report bugs, please e-mail me (Ale) at: pressalex@hotmail.com
