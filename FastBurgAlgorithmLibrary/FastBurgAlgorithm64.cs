﻿using System;

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
        private int _absolutePosition;

        // Naming: 
        // - first letter is the same as in the Koen Vos paper
        // - the part after underscore is my description
        private int _iIterationCounter;
        private int _mCoefficientsNumber;
        private int _nHistoryLengthSamples;
        private readonly double[] _xInputSignal; 
        private double[] _g;
        private double[] _r;
        private double[] _c;
        private double[] _kReflectionCoefs;

        /// <summary>
        /// Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        private double[] _deltaRAndAProduct;

        private double[] _aPredictionCoefs;

        public FastBurgAlgorithm64(double[] inputSignal) 
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
            _aPredictionCoefs = new double[_mCoefficientsNumber + 1];
            _g = new double[_mCoefficientsNumber + 2];
            _r = new double[_mCoefficientsNumber + 1];
            _c = new double[_mCoefficientsNumber + 1];
            _kReflectionCoefs = new double[_mCoefficientsNumber + 1];
            _deltaRAndAProduct = new double[_mCoefficientsNumber + 1];

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
                prediction -= _aPredictionCoefs[index] * 
                    _xInputSignal[_absolutePosition - index];
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
            double[] predictionCoefs = (double[])_aPredictionCoefs.Clone();

            return predictionCoefs;
        }

        /// <summary>
        /// Returns prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double[] GetReflectionCoefs()
        {
            double[] reflectionCoefs = (double[])_kReflectionCoefs.Clone();

            return reflectionCoefs;
        }

        /// <summary>
        /// Updates vector g. For details see step 7 of algorithm on page 3 of 
        /// A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateG()
        {
            double[] oldG = (double[])_g.Clone();

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
                double innerProduct1 = 0;
                double innerProduct2 = 0;
                for (int indexColumn = 0; 
                    indexColumn <= _iIterationCounter; 
                    indexColumn++)
                {
                    innerProduct1 +=
                         _xInputSignal[_absolutePosition - _nHistoryLengthSamples +
                            _iIterationCounter - indexColumn] * 
                        _aPredictionCoefs[indexColumn];
                    innerProduct2 +=
                         _xInputSignal[_absolutePosition - 1 -
                            _iIterationCounter + indexColumn] * 
                        _aPredictionCoefs[indexColumn];
                }

                _deltaRAndAProduct[indexRow] =
                    -_xInputSignal[_absolutePosition - _nHistoryLengthSamples +
                        _iIterationCounter - indexRow] *
                    innerProduct1 -
                    _xInputSignal[_absolutePosition - 1 -
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
            double[] oldR = (double[])_r.Clone();

            for (int index = 0; index <= _iIterationCounter - 1; index++)
            {
                _r[index + 1] = oldR[index] -
                    _xInputSignal[_absolutePosition - _nHistoryLengthSamples + index] *
                    _xInputSignal[_absolutePosition - _nHistoryLengthSamples + _iIterationCounter] -
                    _xInputSignal[_absolutePosition - 1 - index] *
                    _xInputSignal[_absolutePosition - 1 - _iIterationCounter];
            }
            
            _r[0] = 2 * _c[_iIterationCounter + 1];
        }

        /// <summary>
        /// Updates vector of prediction coefficients. For details see step 2 of
        /// algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdatePredictionCoefs()
        {
            double[] oldAPredictionCoefs = (double[])_aPredictionCoefs.Clone();

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
            double nominator = 0;
            double denominator = Double.Epsilon;

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
                Math.Abs(_xInputSignal[_absolutePosition - _nHistoryLengthSamples]) *
                Math.Abs(_xInputSignal[_absolutePosition - _nHistoryLengthSamples]) -
                Math.Abs(_xInputSignal[_absolutePosition - 1]) *
                Math.Abs(_xInputSignal[_absolutePosition - 1]);
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
                    _c[j] += _xInputSignal[index] * _xInputSignal[index + j];
            }
        }
    }
}
