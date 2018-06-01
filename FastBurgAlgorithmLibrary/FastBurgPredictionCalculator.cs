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
            float[] inputAudio,
            float[] forwardPredictions,
            float[] backwardPredictions,
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            Initialization(inputAudio, position, coefficientsNumber, historyLengthSamples);





            double ACCUM = 0.0;

            for (int I = 1; I <= coefficientsNumber; I++)
                ACCUM += inputAudio[position - I] * (-1) * a[I];

            forwardPredictions[position] = (float)ACCUM;

            ACCUM = 0.0;
            for (int I = 1; I <= coefficientsNumber; I++)
                ACCUM += inputAudio[position - historyLengthSamples + I] *
                    (-1) * a[I];

            backwardPredictions[position - historyLengthSamples] =
                (float)ACCUM;
        }

        private static void Initialization(float[] inputAudio, int position, int coefficientsNumber, int historyLengthSamples)
        {
            double[] c = FindAutocorrelation(inputAudio, position, coefficientsNumber, historyLengthSamples);
        }

        private static double[] FindAutocorrelation(float[] inputAudio, int position, int coefficientsNumber, int historyLengthSamples)
        {
            double[] c = new double[coefficientsNumber + 1];

            for (int j = 0; j <= coefficientsNumber; j++)
            {
                c[j] = 0;
                for (int index = position - historyLengthSamples; index <= historyLengthSamples - 1 - j; index++)
                    c[j] += inputAudio[index] * inputAudio[index + j];
            }

            return c;
        }
    }
}
