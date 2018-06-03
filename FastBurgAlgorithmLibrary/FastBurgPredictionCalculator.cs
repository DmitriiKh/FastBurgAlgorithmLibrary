using System;

namespace FastBurgAlgorithmLibrary
{
    public class FastBurgPredictionCalculator
    {
        private int i_iterationCounter;
        private double[,] g;
        private double[] r;

        public double[] a_predictionCoefficients { get; set; }

        public FastBurgPredictionCalculator(int coefficientsNumber)
        {
            a_predictionCoefficients = new double[coefficientsNumber];
            g = new double[coefficientsNumber, 2];
            r = new double[coefficientsNumber];
        }

        /// <summary>
        /// Calculates one prediction error value for one sample using CPU
        /// For details please see "vos_fastburg.pdf"
        /// </summary>
        public void Calculate(
            float[] inputSignal,
            float[] forwardPredictions,
            float[] backwardPredictions,
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            Initialization(inputSignal, position, coefficientsNumber, historyLengthSamples);





           
        }

        private void Initialization(float[] inputSignal, int position, int coefficientsNumber, int historyLengthSamples)
        {
            double[] c = FindAutocorrelation(inputSignal, position, coefficientsNumber, historyLengthSamples);

            i_iterationCounter = 0;
            a_predictionCoefficients[0] = 1;
            g[0, 0] = 2 * c[0] - 
                Math.Pow(Math.Abs(inputSignal[position - historyLengthSamples]), 2) -
                Math.Pow(Math.Abs(inputSignal[position - 1]), 2);
            g[0, 1] = 2 * c[1];
            r[1] = 2 * c[1];
        }

        private double[] FindAutocorrelation(float[] inputSignal, int position, int coefficientsNumber, int historyLengthSamples)
        {
            double[] c = new double[coefficientsNumber + 1];

            for (int j = 0; j <= coefficientsNumber; j++)
            {
                c[j] = 0;
                for (int index = position - historyLengthSamples; index <= historyLengthSamples - 1 - j; index++)
                    c[j] += inputSignal[index] * inputSignal[index + j];
            }

            return c;
        }
    }
}
