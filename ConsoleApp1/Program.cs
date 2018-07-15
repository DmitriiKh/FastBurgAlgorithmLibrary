using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FastBurgAlgorithmLibrary;

namespace ConsoleApp1
{
    class Program
    {
        static void Main(string[] args)
        {
            const int coefNumber = 4;
            const int historyLength = 512;
            const int numberOfSamplesToCheck = 10;

            double[] input_audio =
                new double[historyLength + numberOfSamplesToCheck + 10];

            for (int i = 0; i < input_audio.Length; i++)
            {
                input_audio[i] = System.Math.Sin(
                    2 * System.Math.PI * i / (historyLength / 5.2));
            }

            FastBurgAlgorithm fba = new FastBurgAlgorithm(input_audio, 64);

            for (int index = historyLength;
                index < historyLength + numberOfSamplesToCheck;
                index++)
            {
                fba.Train(index, coefNumber, historyLength);
                var forwardPrediction = fba.GetForwardPrediction();

            }
        }
    }
}
