using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastBurgAlgorithmLibrary
{
    public static class FastBurgPredictionCalculator
    {
        /// <summary>
        /// Calculates one prediction error value for one sample using CPU
        /// For details please see "vos_fastburg.pdf"
        /// </summary>
        public static void Calculate(
            float[] inputaudio,
            float[] forwardPredictions,
            float[] backwardPredictions,
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            Initialization(inputaudio, position, coefficientsNumber, historyLengthSamples);





            double ACCUM = 0.0;

            for (int I = 1; I <= coefficientsNumber; I++)
                ACCUM += inputaudio[position - I] * (-1) * a[I];

            forwardPredictions[position] = (float)ACCUM;

            ACCUM = 0.0;
            for (int I = 1; I <= coefficientsNumber; I++)
                ACCUM += inputaudio[position - historyLengthSamples + I] *
                    (-1) * a[I];

            backwardPredictions[position - historyLengthSamples] =
                (float)ACCUM;
        }

        private static void Initialization(float[] inputaudio, int position, int coefficientsNumber, int historyLengthSamples)
        {
            for (int j = 0; j <= coefficientsNumber; j++)
            {
                for (int )
                c[j] += 5;
            }
        }
    }
}
