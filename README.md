# FastBurgAlgorithmLibrary
Implementation of Fast Burg Algorithm by Koen Vos for real signals (like audio or stocks) in C#.

Master branch unit tests
[![Build Status](https://travis-ci.org/DmitriiKh/FastBurgAlgorithmLibrary.svg?branch=master)]

# Two versions of class
FastBurgAlgorithm128 uses internal variables of type decimal which gives accuracy.
FastBurgAlgorithm64 uses internal variables of type double which gives speed.

# Using:
```csharp
double[] input = new double[2048]; 
    
// for this example we create sinusoid as input signal
for (int i = 0; i < input.Length; i++)
{
    input[i] = System.Math.Sin( 
        2 * System.Math.PI * i / (512 / 5.2));
}

int position = 1234;

// Connect FastBurgAlgorithm instance to input signal
FastBurgAlgorithm128 fba = new FastBurgAlgorithm128(input);
// Train FastBurgAlgorithm at position 1234 with 
// 4 coefficients using 512 previous samples
fba.Train(position, 4, 512);
// get prediction for sample at position 1234
var forwardPrediction = fba.GetForwardPrediction();

```
