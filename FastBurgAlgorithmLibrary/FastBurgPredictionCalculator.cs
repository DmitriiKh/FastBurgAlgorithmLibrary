using System;

namespace FastBurgAlgorithmLibrary
{
    public class FastBurgPredictionCalculator
    {
        private int absolutePosition;
        private int i_iterationCounter;
        private int m_coefficientsNumber;
        private int N_historyLengthSamples;
        private float[] x_inputSignal;
        private double[] g;
        private double[] r;
        private double[] c;
        private double[] k_reflectionCoefs;

        public double[] a_predictionCoefs { get; set; }

        public FastBurgPredictionCalculator(float[] inputSignal)
        {
            x_inputSignal = inputSignal;
            a_predictionCoefs = new double[m_coefficientsNumber];
            g = new double[m_coefficientsNumber + 1];
            r = new double[m_coefficientsNumber];
        }

        /// <summary>
        /// Calculates one prediction error value for one sample using CPU
        /// For details please see "vos_fastburg.pdf"
        /// </summary>
        public void Calculate(
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            Initialization(position, coefficientsNumber, historyLengthSamples);

            ComputeReflectionCoefs();

            UpdatePredictionCoefs();

            i_iterationCounter++;
            if (i_iterationCounter == m_coefficientsNumber)
                return;

            UpdateR();



           
        }

        private void UpdateR()
        {
            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                r[index + 1] = r[index] -
                    x_inputSignal[absolutePosition - N_historyLengthSamples + index] *
                    x_inputSignal[absolutePosition - N_historyLengthSamples + i_iterationCounter] -
                    x_inputSignal[absolutePosition - 1 - index] *
                    x_inputSignal[i_iterationCounter];
            }

            r[0] = 2 * c[i_iterationCounter];
        }

        private void UpdatePredictionCoefs()
        {
            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                a_predictionCoefs[index] = a_predictionCoefs[index] + 
                    k_reflectionCoefs[i_iterationCounter] * 
                    a_predictionCoefs[J_inversOrder(index)];
            }
        }

        private void ComputeReflectionCoefs()
        {
            // for real numbers input signals
            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                k_reflectionCoefs[i_iterationCounter] += g[J_inversOrder(index)] * g[index];
            }
            k_reflectionCoefs[i_iterationCounter] = - k_reflectionCoefs[i_iterationCounter];
        }

        private int J_inversOrder(int index)
        {
            return i_iterationCounter + 1 - index;
        }

        private void Initialization(int position, int coefficientsNumber, int historyLengthSamples)
        {
            m_coefficientsNumber = coefficientsNumber;
            N_historyLengthSamples = historyLengthSamples;
            absolutePosition = position;

            c = FindAutocorrelation();

            i_iterationCounter = 0;
            a_predictionCoefs[0] = 1;
            g[0] = 2 * c[0] - 
                Math.Pow(Math.Abs(x_inputSignal[absolutePosition - historyLengthSamples]), 2) -
                Math.Pow(Math.Abs(x_inputSignal[absolutePosition - 1]), 2);
            g[1] = 2 * c[1];
            // the paper says r[1], error?
            r[0] = 2 * c[1];
        }

        private double[] FindAutocorrelation()
        {
            double[] c = new double[m_coefficientsNumber + 1];

            for (int j = 0; j <= m_coefficientsNumber; j++)
            {
                c[j] = 0;
                for (int index = absolutePosition - N_historyLengthSamples; index <= N_historyLengthSamples - 1 - j; index++)
                    c[j] += x_inputSignal[index] * x_inputSignal[index + j];
            }

            return c;
        }
    }
}
