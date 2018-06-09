using System;

namespace FastBurgAlgorithmLibrary
{
    public class FastBurgPrediction
    {
        private int absolutePosition;
        private int i_iterationCounter;
        private int m_coefficientsNumber;
        private int N_historyLengthSamples;
        private readonly float[] x_inputSignal;
        private double[] g;
        private double[] r;
        private double[] c;
        private double[] k_reflectionCoefs;
        /// <summary>
        /// Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        private double[] deltaRAndAProduct;

        private double[] a_predictionCoefs { get; set; }

        public FastBurgPrediction(float[] inputSignal)
        {
            x_inputSignal = inputSignal;
            a_predictionCoefs = new double[m_coefficientsNumber + 1];
            g = new double[m_coefficientsNumber + 1];
            r = new double[m_coefficientsNumber + 1];
            c = new double[m_coefficientsNumber + 1];
            k_reflectionCoefs = new double[m_coefficientsNumber + 1];
            deltaRAndAProduct = new double[m_coefficientsNumber + 1];
        }

        /// <summary>
        /// Calculates one prediction error value for one sample using CPU
        /// For details please see "vos_fastburg.pdf"
        /// </summary>
        public void Train(
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            absolutePosition = position;

            while (i_iterationCounter < m_coefficientsNumber)
            {
                Initialization(position, coefficientsNumber, historyLengthSamples);

                ComputeReflectionCoefs();

                UpdatePredictionCoefs();

                i_iterationCounter++;
                if (i_iterationCounter == m_coefficientsNumber)
                    return;

                UpdateR();

                ComputeDeltaRMultByA();

                UpdateG();
            }
        }

        public float GetForwardPrediction()
        {
            double prediction = 0;
            for (int index = 0; index <= a_predictionCoefs.Length - 1; index++)
            {
                prediction -= a_predictionCoefs[index] * 
                    x_inputSignal[absolutePosition - 1 - index];
            }

            return (float)prediction;
        }

        private void UpdateG()
        {
            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                g[index] = 
                    g[index] + 
                    k_reflectionCoefs[index] * g[J_inversOrder(index)] + 
                    deltaRAndAProduct[index];
            }

            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                g[i_iterationCounter] += r[index] * a_predictionCoefs[index];
            }
        }

        private void ComputeDeltaRMultByA()
        {
            for (int indexRow = 0; indexRow <= i_iterationCounter; indexRow++)
            {
                double innerProduct1 = 0;
                double innerProduct2 = 0;
                for (int indexColumn = 0; 
                    indexColumn <= i_iterationCounter; 
                    indexColumn++)
                {
                    innerProduct1 += 
                        x_inputSignal[i_iterationCounter - indexColumn] * 
                        a_predictionCoefs[indexColumn];
                    innerProduct2 += 
                        x_inputSignal[absolutePosition - 1 - indexColumn] * 
                        a_predictionCoefs[indexColumn];
                }

                deltaRAndAProduct[indexRow] =
                    -x_inputSignal[indexRow] *
                    innerProduct1 -
                    x_inputSignal[absolutePosition - 1 - indexRow] * 
                    innerProduct2;
            }
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
            // updating r[0] after r[1:i_iterationCounter] are done
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
            double numerator = 0;
            double denominator = 0;
            // for real numbers input signals
            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                numerator += a_predictionCoefs[index] * g[J_inversOrder(index)];
                denominator += a_predictionCoefs[index] * g[index];
            }

            k_reflectionCoefs[i_iterationCounter] = - numerator / denominator;
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
