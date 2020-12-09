import sys
import numpy as np
import math
import matplotlib.pylab as plt

# ====================================
# ===== Generic Parameters ===========
# ====================================
# Observation time interval [seconds]
GP_OBSERVATION_INTERVAL = 0.1

# ====================================
# ===== Digital Filter Parameters ====
# ====================================
# ADC sampling frequency [Hz]
DF_DISCRETIZATION_FREQUENCY = 2000

# Filter type: Low Pass
# Filter model: Butterworth
# Filter order: 2
# Sampling Frequency: 2000 Hz
# Cut Frequency: 200.000000 Hz

# Numerator coefficients
DF_BCOEF = [
    0.06745527606901530200,
    0.13491055213803060000,
    0.06745527606901530200,
       1.4409322,
     -0.88319898,
    0.8232859
]
# De-numerator coefficients
DF_ACOEF = [
    1.00000000000000000000,
    -1.14298050253990090000,
    0.41280159809618860000,
    -0.39832281,
    1.10571257,
    1.1167958
]

# ====================================
# ===== Input Signal Parameters ====+
# ====================================
# Start Frequency [Hz]
INPUT_START_FREQ = 50
# End Frequency [Hz]
INPUT_END_FREQ = 950
# Total Number of Frequencies [Hz]
INPUT_FREQ_POINTS_CNT = 10
# Amplitude of Input Signal [counts]
INPUT_AMPLITUDE = 127

# -------------------------------------------
# Define digital filter order
DF_Order = len(DF_BCOEF) - 1

# Buffer to store digital filter internal data
DF_Buffer = np.zeros(DF_Order + 1)

if (len(DF_ACOEF) < (DF_Order + 1)):
    print('ERROR: DF_BCOEF length is to small!')
    sys.exit()

# This function to process new sample and update
# internal buffer DF_Buffer values
def DigitalFilter_ProcessOneSample(NewSample):
    # Shift the old values in buffer DF_Buffer
    for idx in range(DF_Order, 0, -1):
        DF_Buffer[idx] = DF_Buffer[idx - 1]
    # Sum in Denumerator
    Accumulator = NewSample
    for idx in range(1, DF_Order + 1):
        Accumulator = Accumulator - DF_ACOEF[idx] * DF_Buffer[idx]
    DF_Buffer[0] = Accumulator
    # Sum in Numerator
    Accumulator = 0
    for idx in range(0, DF_Order+1):
        Accumulator = Accumulator + DF_BCOEF[idx] * DF_Buffer[idx]
    return (Accumulator)

# This function to process each sample of input signal by applying
# of the digital filter (call of DigitalFilter_ProcessOneSample)
def ApplyDigitalFilter(SamplesArray):
    total_samples_cnt = len(SamplesArray)
    OutSamples = np.zeros(total_samples_cnt)
    for cur_idx in range(total_samples_cnt):
        OutSamples[cur_idx] = DigitalFilter_ProcessOneSample(SamplesArray[cur_idx])
    return(OutSamples)

# Prepare vector of the input signal frequencies to be tested
input_freq_vector = np.linspace(INPUT_START_FREQ, INPUT_END_FREQ, INPUT_FREQ_POINTS_CNT)
# Time vector according to sampling frequency
time_vector = np.linspace(0, GP_OBSERVATION_INTERVAL, (int)(GP_OBSERVATION_INTERVAL * DF_DISCRETIZATION_FREQUENCY + 1))

k_vect = list()
# Calculate amplitude-frequency response of digital filter by applying of harmonic signals
# to its input with different frequencies
for cur_freq in input_freq_vector:
    # Generate samples of input signal ar current frequency
    input_signal = INPUT_AMPLITUDE * np.sin(2 * math.pi * cur_freq * time_vector)
    # Obtain output samples by applying of the digital filter
    output_signal = ApplyDigitalFilter(input_signal)
    # Remove first sample from output signal to ignore transient process
    output_signal = output_signal[10:-1]
    # Calculate transfer coefficient at current frequency as ratio between output
    # and input amplitudes
    k = (max(output_signal) - min(output_signal)) / (max(input_signal) - min(input_signal))
    # Calculate transfer coefficient in dB and append to list k_vect
    k_vect.append(20 * math.log10(k))
    # Plot figure with output signal at current frequency
    plt.figure()
    plt.plot(output_signal)
    plt.xlabel(f'frequency: {cur_freq} Hz')
# Plot figure with amplitude-frequency response of digital filter
plt.figure()
plt.plot(input_freq_vector, k_vect)
plt.grid()
plt.xlabel(f'Frequency, Hz')
plt.ylabel(f'Transfer coefficient, dB')
plt.show()



