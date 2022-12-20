/* -------------------------------------------------------------------------- *
 *  Libraries
 * -------------------------------------------------------------------------- */

#include <algorithm>
    using std::for_each;

#include <complex>
    using std::complex;
    using std::conj;
    using std::polar;

#include <vector>
    using std::vector;

/* -------------------------------------------------------------------------- *
 *  Declerations
 * -------------------------------------------------------------------------- */

template<typename T>
static void FFT(vector<complex<T>> &sample_);

template<typename T>
static void IFFT(vector<complex<T>> &sample_);

/* -------------------------------------------------------------------------- *
 *  Implementation
 * -------------------------------------------------------------------------- */


template<typename T>
static void FFT(vector<complex<T>> &sample_)
{
    const size_t N = sample_.size() / 2;
    if (N < 1) return;  // Return condition

    // Split samples to even and odd indexes
    vector<complex<T>> even;
    vector<complex<T>> odd;
    for (size_t i = 0; i < N; ++i)
    {
        even.push_back(sample_[2 * i]);
        odd.push_back(sample_[2 * i + 1]);
    }

    // Compute FFT on the divided signals
    FftAlgorithm(even);
    FftAlgorithm(odd);

    // Combine the divided signals back
    for (size_t k = 0; k < N; ++k)
    {
        complex<T> odd_factor = polar(T(1), -T(M_PI) * k / N) * odd[k];
        sample_[k] = even[k] + odd_factor;
        sample_[k + N] = even[k] - odd_factor;
    }
}

template<typename T>
static void IFFT(vector<complex<T>> &sample_)
{
    const size_t N = sample_.size();

    for_each(sample_.begin(), sample_.end(),
        [](complex<T> &val_){ val_ = conj(val_); });

    FftAlgorithm(sample_);

    for_each(sample_.begin(), sample_.end(),
        [N](complex<T> &val_){ val_ = conj(val_) / T(N); });

}