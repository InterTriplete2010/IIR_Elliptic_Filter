#pragma once

#include <stdio.h> 
#include <complex.h>

//#define ARMA_DONT_USE_CXX11
#define ARMA_DONT_USE_CXX11_MUTEX
#include <armadillo>

#ifndef Elliptic_Filter_H
#define Elliptic_Filter_H

#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif

    namespace IIR_E_F
    {

        class IIR_Elliptic_Filter

        {


        private:

            //get analog, pre - warped frequencies
            void freq_pre_wrapped(int, double, double);

            //convert to low-pass prototype estimate
            void Wn_f1_Wn_f2(int, double, double);

            //Get N - th order Elliptic analog lowpass prototype. 
            void ellipap(int, double, double);

            //Get Landen vector of descending moduli and return a vector of double numbers
            std::vector<double> landen_vector(double);

            //Get Landen vector of descending moduli and return a vector of complex numbers
            std::vector<std::complex<double>> landen_vector_complex(double);

            //Elliptic function with normalized complex argument
            double sne(std::vector<double>, double);

            //Elliptic function with normalized complex argument for complex values
            std::complex<double> sne(std::complex<double>, double);

            //Elliptic function with normalized complex argument
            std::vector<double> cde(std::vector<double>, double);

            //Elliptic function with normalized complex argument for complex values
            std::vector<std::complex<double>> cde(std::vector<std::complex<double>>, double);

            //Inverse of sn elliptic function
            std::complex<double> asne(std::complex<double>, double);

            //Complete elliptic integral of first kind
            std::vector<double> ellipk(double);

            //Sort complex numbers into complex conjugate pairs, starting from with the number with the lowest real part.
            //If the real part is the same for all the number, order according to the absolute value of the highest imaginary part.
            //Within a pair, the element with negative imaginary part comes first. 
            std::vector<std::complex<double>> cplxpair(std::vector<std::complex<double>>);

            //Transform to state-space
            void zp2ss();
           
            //Bilinear transformation to find discrete equivalent
            void bilinear(arma::mat, arma::mat, arma::mat, arma::mat, double, int, int);

            //Calculate the coefficients of the polynomial (based on Matlab code)
            std::vector<std::complex<double>> poly(std::vector<std::complex<double>>, int);

            //Calculate the coefficients of the polynomial (based on Matlab code)
            std::vector<double> poly(arma::mat, int);

            //Calculate the coefficients of the polynomial (based on Matlab code)
            std::vector<double> poly(arma::cx_vec, int);

            //Extract the zeros of the state-space system
            void sss_zeros();

        public:

            //Estimate the coeffients of a band-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2bp(int, double, double, double, double);

            //Estimate the coeffients of a band-stop filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2bs(int, double, double, double, double);

            //Estimate the coeffients of a high-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2hp(int, double, double, double);

            //Estimate the coeffients of a low-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2lp(int, double, double, double);

            //Check the stability of the filter. Returns "true" is the filter is stable, false if it is unstable 
            bool check_stability_iir(std::vector<std::vector<double> >);

            //Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
            std::vector<double> Filter_Data(std::vector<std::vector<double> > coeff_filt, std::vector<double> pre_filt_signal);

        };

    }

#endif

#ifdef __cplusplus

}

#endif

