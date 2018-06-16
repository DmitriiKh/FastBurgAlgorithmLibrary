using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Fast implimentation of Burg algorithm for real signals.
    /// For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    /// FastBurgAlgorithm64 uses internal variables of type double
    /// </summary>
    public class FastBurgAlgorithm64
    {
        /// <summary>
        /// Position in x_inputSignal that we need prediction for. 
        /// </summary>
        private int absolutePosition;

        // Naming: 
        // - first letter is the same as in the Koen Vos paper
        // - the part after underscore is my description
        private int i_iterationCounter;
        private int m_coefficientsNumber;
        private int N_historyLengthSamples;
        private readonly double[] x_inputSignal; 
        private double[] g;
        private double[] r;
        private double[] c;
        private double[] k_reflectionCoefs;

        /// <summary>
        /// Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        private double[] deltaRAndAProduct;

        private double[] a_predictionCoefs;

        public FastBurgAlgorithm64(double[] inputSignal) 
        {
            x_inputSignal = inputSignal;
        }

        /// <summary>
        /// Calculates prediction coefficients for one sample using CPU
        /// </summary>
        /// <param name="position"> Position in inputSignal that we need 
        /// prediction for. Must be greater than historyLengthSamples </param>
        /// <param name="coefficientsNumber"> Number of prediction coefficients
        /// that will be calculated. Greater number gives more accurate 
        /// prediction but takes more time to calculate </param>
        /// <param name="historyLengthSamples"> Number of samples that will 
        /// be used to calculate prediction coefficients </param>
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

        /// <summary>
        /// Returns forward prediction based on prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double GetForwardPrediction() 
        {
            double prediction = 0;
            for (int index = 1; index <= a_predictionCoefs.Length - 1; index++)
            {
                prediction -= a_predictionCoefs[index] * 
                    x_inputSignal[absolutePosition - index];
            }

            return prediction; 
        }

        /// <summary>
        /// Returns prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double[] GetPredictionCoefs()
        {
            double[] predictionCoefs = (double[])a_predictionCoefs.Clone();

            return predictionCoefs;
        }

        /// <summary>
        /// Returns prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double[] GetReflectionCoefs()
        {
            double[] reflectionCoefs = (double[])k_reflectionCoefs.Clone();

            return reflectionCoefs;
        }

        /// <summary>
        /// Updates vector g. For details see step 7 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateG()
        {
            double[] old_g = (double[])g.Clone();

            // g.Length is i_iterationCounter + 1
            for (int index = 0; index <= i_iterationCounter; index++)
            {
                g[index] = 
                    old_g[index] + 
                    k_reflectionCoefs[i_iterationCounter - 1] * old_g[J_inversOrder(index, i_iterationCounter)] + 
                    deltaRAndAProduct[index];
            }

            for (int index = 0; index <= i_iterationCounter; index++)
            {
                g[i_iterationCounter + 1] += r[index] * a_predictionCoefs[index];
            }
        }

        /// <summary>
        /// Calculates vector deltaRAndAProduct. For details see step 6 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
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

        /// <summary>
        /// Updates vector r. For details see step 5 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateR()
        {
            double[] old_r = (double[])r.Clone();

            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                r[index + 1] = old_r[index] -
                    x_inputSignal[absolutePosition - N_historyLengthSamples + index] *
                    x_inputSignal[absolutePosition - N_historyLengthSamples + i_iterationCounter] -
                    x_inputSignal[absolutePosition - 1 - index] *
                    x_inputSignal[absolutePosition - 1 - i_iterationCounter];
            }
            
            r[0] = 2 * c[i_iterationCounter + 1];
        }

        /// <summary>
        /// Updates vector of prediction coefficients. For details see step 2 of
        /// algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdatePredictionCoefs()
        {
            double[] old_a_predictionCoefs = (double[])a_predictionCoefs.Clone();

            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                a_predictionCoefs[index] = old_a_predictionCoefs[index] + 
                    k_reflectionCoefs[i_iterationCounter] * 
                    old_a_predictionCoefs[J_inversOrder(index, i_iterationCounter + 1)];
            }
        }

        /// <summary>
        /// Computes vector of reflection coefficients. For details see step 1
        /// of algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
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
        /// Inverts index to flip a vector insted of multiplication with J matrix.
        /// For details see (12) on page 2 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        /// <param name="index">from 0 to max</param>
        /// <param name="max">positive number</param>
        /// <returns></returns>
        private int J_inversOrder(int index, int max)
        {
            return max - index;
        }

        /// <summary>
        /// Initializes i_iterationCounter and vectors. For details see step 0 of
        /// algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void Initialization()
        {
            FindAutocorrelations();

            i_iterationCounter = 0;
            a_predictionCoefs[0] = 1;
            g[0] = 2 * c[0] -
                Math.Abs(x_inputSignal[absolutePosition - N_historyLengthSamples]) *
                Math.Abs(x_inputSignal[absolutePosition - N_historyLengthSamples]) -
                Math.Abs(x_inputSignal[absolutePosition - 1]) *
                Math.Abs(x_inputSignal[absolutePosition - 1]);
            g[1] = 2 * c[1];
            // the paper says r[1], error in paper?
            r[0] = 2 * c[1];
        }

        /// <summary>
        /// Calculates autocorrelations. For details see step 0 of
        /// algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
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
