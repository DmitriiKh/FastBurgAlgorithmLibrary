using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Implimentation of fast Burg algorithm
    /// </summary>
    public class FastBurgAlgorithm
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

        private double[] a_predictionCoefs;

        public FastBurgAlgorithm(float[] inputSignal)
        {
            x_inputSignal = inputSignal;
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
            m_coefficientsNumber = coefficientsNumber;
            N_historyLengthSamples = historyLengthSamples;
            a_predictionCoefs = new double[m_coefficientsNumber + 1];
            g = new double[m_coefficientsNumber + 2];
            r = new double[m_coefficientsNumber + 1];
            c = new double[m_coefficientsNumber + 1];
            k_reflectionCoefs = new double[m_coefficientsNumber + 1];
            deltaRAndAProduct = new double[m_coefficientsNumber + 1];

            Initialization();

            while (i_iterationCounter <= m_coefficientsNumber)
            {
                ComputeReflectionCoef();

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
            for (int index = 1; index <= a_predictionCoefs.Length - 1; index++)
            {
                prediction -= a_predictionCoefs[index] * 
                    x_inputSignal[absolutePosition - index];
            }

            return (float)prediction;
        }

        private void UpdateG()
        {
            double[] old_g = new double[g.Length];
            for (int index = 0; index < g.Length; index++)
                old_g[index] = g[index];

            // g.Length is i_iterationCounter + 1
            for (int index = 0; index <= i_iterationCounter /*+ 1*/; index++)
            {
                g[index] = 
                    old_g[index] + 
                    k_reflectionCoefs[i_iterationCounter - 1/*index*/] * old_g[J_inversOrder(index, i_iterationCounter/* + 1*/)] + 
                    deltaRAndAProduct[index];
            }

            for (int index = 0; index <= i_iterationCounter; index++)
            {
                g[i_iterationCounter + 1/*2*/] += r[index] * a_predictionCoefs[index];
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
                        x_inputSignal[absolutePosition - N_historyLengthSamples + 
                            i_iterationCounter - indexColumn] * 
                        a_predictionCoefs[indexColumn];
                    innerProduct2 += 
                        x_inputSignal[absolutePosition - 1 - 
                            i_iterationCounter + indexColumn] * 
                        a_predictionCoefs[indexColumn];
                }

                deltaRAndAProduct[indexRow] =
                    -x_inputSignal[absolutePosition - N_historyLengthSamples + 
                        i_iterationCounter - indexRow] *
                    innerProduct1 -
                    x_inputSignal[absolutePosition - 1 - 
                        i_iterationCounter + indexRow] * 
                    innerProduct2;
            }
        }

        private void UpdateR()
        {
            double[] old_r = new double[r.Length];
            for (int index = 0; index < r.Length; index++)
                old_r[index] = r[index];

            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                r[index + 1] = old_r[index] -
                    x_inputSignal[absolutePosition - N_historyLengthSamples + index] *
                    x_inputSignal[absolutePosition - N_historyLengthSamples + i_iterationCounter] -
                    x_inputSignal[absolutePosition - 1 - index] *
                    x_inputSignal[absolutePosition - 1 - i_iterationCounter];
            }
            // updating r[0] after r[1:i_iterationCounter] are done
            r[0] = 2 * c[i_iterationCounter + 1];
        }

        private void UpdatePredictionCoefs()
        {
            double[] old_a_predictionCoefs = new double[a_predictionCoefs.Length];
            for (int index = 0; index < a_predictionCoefs.Length; index++)
                old_a_predictionCoefs[index] = a_predictionCoefs[index];

            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                a_predictionCoefs[index] = old_a_predictionCoefs[index] + 
                    k_reflectionCoefs[i_iterationCounter] * 
                    old_a_predictionCoefs[J_inversOrder(index, i_iterationCounter + 1)];
            }
        }

        private void ComputeReflectionCoef()
        {
            double nominator = 0;
            double denominator = 0;

            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                nominator += a_predictionCoefs[index] * 
                    g[J_inversOrder(index, i_iterationCounter + 1)];
                denominator += a_predictionCoefs[index] * g[index];
            }

            k_reflectionCoefs[i_iterationCounter] = - nominator / denominator;
        }

        /// <summary>
        /// Inverts index to flip a vector 
        /// </summary>
        /// <param name="index">from 0 to max</param>
        /// <param name="max">from 0</param>
        /// <returns></returns>
        private int J_inversOrder(int index, int max)
        {
            return max - index;
        }

        private void Initialization()
        {
            FindAutocorrelations();

            i_iterationCounter = 0;
            a_predictionCoefs[0] = 1;
            g[0] = 2 * c[0] - 
                Math.Pow(Math.Abs(x_inputSignal[absolutePosition - N_historyLengthSamples]), 2) -
                Math.Pow(Math.Abs(x_inputSignal[absolutePosition - 1]), 2);
            g[1] = 2 * c[1];
            // the paper says r[1], error?
            r[0] = 2 * c[1];
        }

        private void FindAutocorrelations()
        {
            for (int j = 0; j <= m_coefficientsNumber; j++)
            {
                c[j] = 0;
                for (int index = absolutePosition - N_historyLengthSamples; 
                    index <= absolutePosition - 1 - j; 
                    index++)
                    c[j] += x_inputSignal[index] * x_inputSignal[index + j];
            }
        }
    }
}
