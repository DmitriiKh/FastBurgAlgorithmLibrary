using System;

namespace FastBurgAlgorithmLibrary
{
    public class FastBurgPredictionCalculator
    {
        private int i_iterationCounter;
        private double[] g;
        private double[] r;
        private double[] k_reflectionCoefs;

        public double[] a_predictionCoefs { get; set; }

        public FastBurgPredictionCalculator(int m_coefficientsNumber)
        {
            a_predictionCoefs = new double[m_coefficientsNumber];
            g = new double[m_coefficientsNumber + 1];
            r = new double[m_coefficientsNumber];
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

            ComputeReflectionCoefs();

            UpdatePredictionCoefs();



           
        }

        private void UpdatePredictionCoefs()
        {
            throw new NotImplementedException();
        }

        private void ComputeReflectionCoefs()
        {
            // for real numbers input signals
            for (int index = 0; index <= i_iterationCounter + 1)
            {
                k_reflectionCoefs[i_iterationCounter] += g[J_inversOrder(index)] * g[index];
            }
            k_reflectionCoefs[i_iterationCounter] = - k_reflectionCoefs[i_iterationCounter];
        }

        private int J_inversOrder(int index)
        {
            return i_iterationCounter + 1 - index;
        }

        private void Initialization(float[] inputSignal, int position, int coefficientsNumber, int historyLengthSamples)
        {
            double[] c = FindAutocorrelation(inputSignal, position, coefficientsNumber, historyLengthSamples);

            i_iterationCounter = 0;
            a_predictionCoefs[0] = 1;
            g[0] = 2 * c[0] - 
                Math.Pow(Math.Abs(inputSignal[position - historyLengthSamples]), 2) -
                Math.Pow(Math.Abs(inputSignal[position - 1]), 2);
            g[1] = 2 * c[1];
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
