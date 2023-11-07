# FastBurgAlgorithmLibrary
This is an implementation of the Fast Burg Algorithm by Koen Vos for real signals (like audio or stocks) in C#.

![Build Status](https://github.com/DmitriiKh/FastBurgAlgorithmLibrary/actions/workflows/dotnet.yml/badge.svg)

# Two versions of the class
FastBurgAlgorithm128 uses internal variables of type decimal which gives more accuracy.
FastBurgAlgorithm64 uses internal variables of type double which gives more speed.

# Using:
```csharp
// The position in which we are looking for a prediction (position > window)
int position = 1234;
// The number of prediction coefficients we are using (larger is more accurate but slower)
int coefsNumber = 4;
// The number of samples (always 2^n) which will be used for producing prediction coefficients
int window = 512;
// The input signal array must contain at least 512 (window) samples before the position and 4 (coefsNumber) after it
// [...optional...] + [window] + [position] + [coefsNumber] + [...optional...]
double[] input = new double[2048]; 
    
// For this example we create a sinusoid as the input signal
for (int i = 0; i < input.Length; i++)
{
    input[i] = System.Math.Sin( 
        2 * System.Math.PI * i / (512 / 5.2));
}

// Pass the input signal to a FastBurgAlgorithm instance 
FastBurgAlgorithm64 fba = new FastBurgAlgorithm64(input);
// Train FastBurgAlgorithm at position 1234 using 512 previous samples and produce 4 prediction coefficients 
fba.Train(position, coefsNumber, window);
// [Optional] Get prediction coefficients
double[] predictionCoefs = fba.GetPredictionCoefs();
// Get a prediction for a sample at position 1234
double forwardPrediction = fba.GetForwardPrediction();

```
