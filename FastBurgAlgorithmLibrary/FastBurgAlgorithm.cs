using System;
using System.Linq;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Fast implimentation of Burg algorithm for real signals.
    /// For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    /// FastBurgAlgorithm uses internal variables of dynamic type 
    /// </summary>
    public abstract class FastBurgAlgorithm
    {
        /// <summary>
        /// Position in x_inputSignal that we need prediction for. 
        /// </summary>
        private int absolutePosition;

        internal readonly Type typeOfInternalVariables;

        // Naming: 
        // - first letter is the same as in the Koen Vos paper
        // - the part after underscore is my description
        private int i_iterationCounter;
        private int m_coefficientsNumber;
        private int N_historyLengthSamples;
        private readonly double[] x_inputSignal; 
        private dynamic g;
        private dynamic r;
        private dynamic c;
        internal dynamic k_reflectionCoefs;
        internal dynamic a_predictionCoefs;

        /// <summary>
        /// Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        private dynamic deltaRAndAProduct;

        protected FastBurgAlgorithm(double[] inputSignal, Type type) 
        {
            x_inputSignal = inputSignal;
            typeOfInternalVariables = type;
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

            CreateInternalVariables();

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
        /// Creates internal variables
        /// </summary>
        private void CreateInternalVariables()
        {
            int length = m_coefficientsNumber + 1;
            a_predictionCoefs =
                Array.CreateInstance(typeOfInternalVariables, length + 1);
            g = Array.CreateInstance(typeOfInternalVariables, length + 1);
            r = Array.CreateInstance(typeOfInternalVariables, length);
            c = Array.CreateInstance(typeOfInternalVariables, length);
            k_reflectionCoefs = 
                Array.CreateInstance(typeOfInternalVariables, length);
            deltaRAndAProduct = 
                Array.CreateInstance(typeOfInternalVariables, length);
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
                prediction -= (double) a_predictionCoefs[index] * 
                    x_inputSignal[absolutePosition - index];
            }

            return prediction; 
        }
        
        /// <summary>
        /// Updates vector g. For details see step 7 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateG()
        {
            decimal[] old_g = (decimal[])g.Clone();

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
                decimal innerProduct1 = 0;
                decimal innerProduct2 = 0;
                for (int indexColumn = 0; 
                    indexColumn <= i_iterationCounter; 
                    indexColumn++)
                {
                    innerProduct1 += 
                        (decimal) x_inputSignal[absolutePosition - N_historyLengthSamples + 
                            i_iterationCounter - indexColumn] * 
                        a_predictionCoefs[indexColumn];
                    innerProduct2 += 
                        (decimal) x_inputSignal[absolutePosition - 1 - 
                            i_iterationCounter + indexColumn] * 
                        a_predictionCoefs[indexColumn];
                }

                deltaRAndAProduct[indexRow] =
                    -(decimal)x_inputSignal[absolutePosition - N_historyLengthSamples + 
                        i_iterationCounter - indexRow] *
                    innerProduct1 -
                    (decimal)x_inputSignal[absolutePosition - 1 - 
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
            decimal[] old_r = (decimal[])r.Clone();

            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                r[index + 1] = old_r[index] -
                    (decimal)x_inputSignal[absolutePosition - N_historyLengthSamples + index] *
                    (decimal)x_inputSignal[absolutePosition - N_historyLengthSamples + i_iterationCounter] -
                    (decimal)x_inputSignal[absolutePosition - 1 - index] *
                    (decimal)x_inputSignal[absolutePosition - 1 - i_iterationCounter];
            }
            
            r[0] = 2 * c[i_iterationCounter + 1];
        }

        /// <summary>
        /// Updates vector of prediction coefficients. For details see step 2 of
        /// algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdatePredictionCoefs()
        {
            decimal[] old_a_predictionCoefs = (decimal[])a_predictionCoefs.Clone();

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
            decimal nominator = 0;
            decimal denominator = 0;

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
                Math.Abs((decimal)x_inputSignal[absolutePosition - N_historyLengthSamples]) *
                Math.Abs((decimal)x_inputSignal[absolutePosition - N_historyLengthSamples]) -
                Math.Abs((decimal)x_inputSignal[absolutePosition - 1]) *
                Math.Abs((decimal)x_inputSignal[absolutePosition - 1]);
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
