using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Fast implimentation of Burg algorithm for real signals.
    /// For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    /// FastBurgAlgorithm128 uses internal variables of type decimal
    /// </summary>
    public class FastBurgAlgorithm128
    {
        /// <summary>
        /// Position in x_inputSignal that we need prediction for. 
        /// </summary>
        private int _absolutePosition;

        // Naming: 
        // - first letter is the same as in the Koen Vos paper
        // - the part after is my description
        private int _iIterationCounter;
        private int _mCoefficientsNumber;
        private int _nHistoryLengthSamples;
        private readonly double[] _xInputSignal; 
        private decimal[] _g;
        private decimal[] _r;
        private decimal[] _c;
        private decimal[] _kReflectionCoefs;

        /// <summary>
        /// Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        private decimal[] _deltaRAndAProduct;

        private decimal[] _aPredictionCoefs;

        public FastBurgAlgorithm128(double[] inputSignal) 
        {
            _xInputSignal = inputSignal;
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
            _absolutePosition = position;
            _mCoefficientsNumber = coefficientsNumber;
            _nHistoryLengthSamples = historyLengthSamples;
            _aPredictionCoefs = new decimal[_mCoefficientsNumber + 1];
            _g = new decimal[_mCoefficientsNumber + 2];
            _r = new decimal[_mCoefficientsNumber + 1];
            _c = new decimal[_mCoefficientsNumber + 1];
            _kReflectionCoefs = new decimal[_mCoefficientsNumber + 1];
            _deltaRAndAProduct = new decimal[_mCoefficientsNumber + 1];

            Initialization();

            while (_iIterationCounter <= _mCoefficientsNumber)
            {
                ComputeReflectionCoef();

                UpdatePredictionCoefs();

                _iIterationCounter++;
                if (_iIterationCounter == _mCoefficientsNumber)
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
            for (int index = 1; index <= _aPredictionCoefs.Length - 1; index++)
            {
                prediction -= (double) _aPredictionCoefs[index] * 
                    _xInputSignal[_absolutePosition - index];
            }

            return prediction; 
        }

        /// <summary>
        /// Returns prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public decimal[] GetPredictionCoefs()
        {
            decimal[] predictionCoefs = (decimal[])_aPredictionCoefs.Clone();

            return predictionCoefs;
        }

        /// <summary>
        /// Returns prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public decimal[] GetReflectionCoefs()
        {
            decimal[] reflectionCoefs = (decimal[])_kReflectionCoefs.Clone();

            return reflectionCoefs;
        }

        /// <summary>
        /// Updates vector g. For details see step 7 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateG()
        {
            decimal[] oldG = (decimal[])_g.Clone();

            // g.Length is i_iterationCounter + 1
            for (int index = 0; index <= _iIterationCounter; index++)
            {
                _g[index] = 
                    oldG[index] + 
                    _kReflectionCoefs[_iIterationCounter - 1] * oldG[J_inversOrder(index, _iIterationCounter)] + 
                    _deltaRAndAProduct[index];
            }

            for (int index = 0; index <= _iIterationCounter; index++)
            {
                _g[_iIterationCounter + 1] += _r[index] * _aPredictionCoefs[index];
            }
        }

        /// <summary>
        /// Calculates vector deltaRAndAProduct. For details see step 6 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void ComputeDeltaRMultByA()
        {
            for (int indexRow = 0; indexRow <= _iIterationCounter; indexRow++)
            {
                decimal innerProduct1 = 0;
                decimal innerProduct2 = 0;
                for (int indexColumn = 0; 
                    indexColumn <= _iIterationCounter; 
                    indexColumn++)
                {
                    innerProduct1 += 
                        (decimal) _xInputSignal[_absolutePosition - _nHistoryLengthSamples + 
                            _iIterationCounter - indexColumn] * 
                        _aPredictionCoefs[indexColumn];
                    innerProduct2 += 
                        (decimal) _xInputSignal[_absolutePosition - 1 - 
                            _iIterationCounter + indexColumn] * 
                        _aPredictionCoefs[indexColumn];
                }

                _deltaRAndAProduct[indexRow] =
                    -(decimal)_xInputSignal[_absolutePosition - _nHistoryLengthSamples + 
                        _iIterationCounter - indexRow] *
                    innerProduct1 -
                    (decimal)_xInputSignal[_absolutePosition - 1 - 
                        _iIterationCounter + indexRow] * 
                    innerProduct2;
            }
        }

        /// <summary>
        /// Updates vector r. For details see step 5 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateR()
        {
            decimal[] oldR = (decimal[])_r.Clone();

            for (int index = 0; index <= _iIterationCounter - 1; index++)
            {
                _r[index + 1] = oldR[index] -
                    (decimal)_xInputSignal[_absolutePosition - _nHistoryLengthSamples + index] *
                    (decimal)_xInputSignal[_absolutePosition - _nHistoryLengthSamples + _iIterationCounter] -
                    (decimal)_xInputSignal[_absolutePosition - 1 - index] *
                    (decimal)_xInputSignal[_absolutePosition - 1 - _iIterationCounter];
            }
            
            _r[0] = 2 * _c[_iIterationCounter + 1];
        }

        /// <summary>
        /// Updates vector of prediction coefficients. For details see step 2 of
        /// algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdatePredictionCoefs()
        {
            decimal[] oldAPredictionCoefs = (decimal[])_aPredictionCoefs.Clone();

            for (int index = 0; index <= _iIterationCounter + 1; index++)
            {
                _aPredictionCoefs[index] = oldAPredictionCoefs[index] + 
                    _kReflectionCoefs[_iIterationCounter] * 
                    oldAPredictionCoefs[J_inversOrder(index, _iIterationCounter + 1)];
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

            for (int index = 0; index <= _iIterationCounter + 1; index++)
            {
                nominator += _aPredictionCoefs[index] * 
                    _g[J_inversOrder(index, _iIterationCounter + 1)];
                denominator += _aPredictionCoefs[index] * _g[index];
            }

            _kReflectionCoefs[_iIterationCounter] = - nominator / denominator;
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

            _iIterationCounter = 0;
            _aPredictionCoefs[0] = 1;
            _g[0] = 2 * _c[0] -
                Math.Abs((decimal)_xInputSignal[_absolutePosition - _nHistoryLengthSamples]) *
                Math.Abs((decimal)_xInputSignal[_absolutePosition - _nHistoryLengthSamples]) -
                Math.Abs((decimal)_xInputSignal[_absolutePosition - 1]) *
                Math.Abs((decimal)_xInputSignal[_absolutePosition - 1]);
            _g[1] = 2 * _c[1];
            // the paper says r[1], error in paper?
            _r[0] = 2 * _c[1];
        }

        /// <summary>
        /// Calculates autocorrelations. For details see step 0 of
        /// algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void FindAutocorrelations()
        {
            for (int j = 0; j <= _mCoefficientsNumber; j++)
            {
                _c[j] = 0;
                for (int index = _absolutePosition - _nHistoryLengthSamples; 
                    index <= _absolutePosition - 1 - j; 
                    index++)
                    _c[j] += (decimal) _xInputSignal[index] * (decimal) _xInputSignal[index + j];
            }
        }
    }
}
