# FastBurgAlgorithmLibrary
Implementation of Fast Burg Algorithm by Koen Vos for real signals (like audio or stocks) in C#.

Master branch unit tests
[![Build Status](https://travis-ci.org/DmitriiKh/FastBurgAlgorithmLibrary.svg?branch=master)]


# Using:
```
double[] input_audio = new double[2048]; 
    
// for this example we create sinusoid as input signal
for (int i = 0; i < input_audio.Length; i++)
{
    input_audio[i] = System.Math.Sin( 
        2 * System.Math.PI * i / (historyLength / 5.2));
}

int position = 1234;

// Connect FastBurgAlgorithm instance to input signal
FastBurgAlgorithm fba = new FastBurgAlgorithm(input_audio);
// Train FastBurgAlgorithm at position 1234 with 
// 4 coefficients using 512 previous samples
fba.Train(position, 4, 512);
// get prediction for sample at position 1234
var forwardPrediction = fba.GetForwardPrediction();

```
