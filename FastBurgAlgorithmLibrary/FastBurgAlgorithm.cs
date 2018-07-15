using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Fast implimentation of Burg algorithm for real signals.
    /// For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    /// FastBurgAlgorithm uses internal variables of dynamic type 
    /// </summary>
    public class FastBurgAlgorithm
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

        private Data _data;

        public FastBurgAlgorithm(double[] inputSignal) 
        {
            x_inputSignal = inputSignal;
        }

        /// <summary>
        /// Calculates prediction coefficients for one sample using CPU
        /// </summary>
        /// <param name="precision">64 or 128, sets internal variables length</param>
        /// <param name="position"> Position in inputSignal that we need 
        /// prediction for. Must be greater than historyLengthSamples </param>
        /// <param name="coefficientsNumber"> Number of prediction coefficients
        /// that will be calculated. Greater number gives more accurate 
        /// prediction but takes more time to calculate </param>
        /// <param name="historyLengthSamples"> Number of samples that will 
        /// be used to calculate prediction coefficients </param>
        public void Train(
            int precision,
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            absolutePosition = position;
            m_coefficientsNumber = coefficientsNumber;
            N_historyLengthSamples = historyLengthSamples;

            CreateInternalVariables(precision);

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

        private void CreateInternalVariables(int precision)
        {
            switch (precision)
            {
                case 64:
                    _data = new Data64(m_coefficientsNumber);
                    break;
                case 128:
                    _data = new Data128(m_coefficientsNumber);
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
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
            for (int index = 1; index <= _data.a_predictionCoefs.Length - 1; index++)
            {
                prediction -= (double) _data.a_predictionCoefs[index] * 
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
            Array.Copy(_data.g, _data.old_g, m_coefficientsNumber + 2);

            // g.Length is i_iterationCounter + 1
            for (int index = 0; index <= i_iterationCounter; index++)
            {
                _data.g[index] = 
                    _data.old_g[index] +
                    _data.k_reflectionCoefs[i_iterationCounter - 1] * _data.old_g[J_inversOrder(index, i_iterationCounter)] +
                    _data.deltaRAndAProduct[index];
            }

            for (int index = 0; index <= i_iterationCounter; index++)
            {
                _data.g[i_iterationCounter + 1] += _data.r[index] * _data.a_predictionCoefs[index];
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
                dynamic innerProduct1 = 0;
                dynamic innerProduct2 = 0;
                for (int indexColumn = 0; 
                    indexColumn <= i_iterationCounter; 
                    indexColumn++)
                {
                    innerProduct1 +=
                        _data.DoubleToInternal(x_inputSignal[absolutePosition - N_historyLengthSamples + 
                            i_iterationCounter - indexColumn]) * 
                        _data.a_predictionCoefs[indexColumn];
                    innerProduct2 +=
                        _data.DoubleToInternal(x_inputSignal[absolutePosition - 1 - 
                            i_iterationCounter + indexColumn]) * 
                        _data.a_predictionCoefs[indexColumn];
                }

                _data.deltaRAndAProduct[indexRow] =
                    -_data.DoubleToInternal(x_inputSignal[absolutePosition - N_historyLengthSamples + 
                        i_iterationCounter - indexRow]) *
                    innerProduct1 -
                    _data.DoubleToInternal(x_inputSignal[absolutePosition - 1 - 
                        i_iterationCounter + indexRow]) * 
                    innerProduct2;
            }
        }

        /// <summary>
        /// Updates vector r. For details see step 5 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateR()
        {
            Array.Copy(_data.r, _data.old_r, m_coefficientsNumber + 1);

            for (int index = 0; index <= i_iterationCounter - 1; index++)
            {
                _data.r[index + 1] = _data.old_r[index] -
                    _data.DoubleToInternal(x_inputSignal[absolutePosition - N_historyLengthSamples + index]) *
                    _data.DoubleToInternal(x_inputSignal[absolutePosition - N_historyLengthSamples + i_iterationCounter]) -
                                     _data.DoubleToInternal(x_inputSignal[absolutePosition - 1 - index]) *
                                     _data.DoubleToInternal(x_inputSignal[absolutePosition - 1 - i_iterationCounter]);
            }

            _data.r[0] = 2 * _data.c[i_iterationCounter + 1];
        }

        /// <summary>
        /// Updates vector of prediction coefficients. For details see step 2 of
        /// algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdatePredictionCoefs()
        {
            Array.Copy(_data.a_predictionCoefs, _data.old_a_predictionCoefs, m_coefficientsNumber + 1);

            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                _data.a_predictionCoefs[index] = _data.old_a_predictionCoefs[index] +
                                                _data.k_reflectionCoefs[i_iterationCounter] * 
                    _data.old_a_predictionCoefs[J_inversOrder(index, i_iterationCounter + 1)];
            }
        }

        /// <summary>
        /// Computes vector of reflection coefficients. For details see step 1
        /// of algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void ComputeReflectionCoef()
        {
            dynamic nominator = 0;
            dynamic denominator = 0;

            for (int index = 0; index <= i_iterationCounter + 1; index++)
            {
                nominator += _data.a_predictionCoefs[index] *
                             _data.g[J_inversOrder(index, i_iterationCounter + 1)];
                denominator += _data.a_predictionCoefs[index] * _data.g[index];
            }

            _data.k_reflectionCoefs[i_iterationCounter] = - nominator / denominator;
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
            _data.a_predictionCoefs[0] = 1;
            _data.g[0] = 2 * _data.c[0] -
                Math.Abs(_data.DoubleToInternal(x_inputSignal[absolutePosition - N_historyLengthSamples])) *
                Math.Abs(_data.DoubleToInternal(x_inputSignal[absolutePosition - N_historyLengthSamples])) -
                Math.Abs(_data.DoubleToInternal(x_inputSignal[absolutePosition - 1])) *
                Math.Abs(_data.DoubleToInternal(x_inputSignal[absolutePosition - 1]));
            _data.g[1] = 2 * _data.c[1];
            // the paper says r[1], error in paper?
            _data.r[0] = 2 * _data.c[1];
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
                _data.c[j] = 0;
                for (int index = absolutePosition - N_historyLengthSamples; 
                    index <= absolutePosition - 1 - j; 
                    index++)
                    _data.c[j] += _data.DoubleToInternal(x_inputSignal[index]) * _data.DoubleToInternal(x_inputSignal[index + j]);
            }
        }
    }
}
