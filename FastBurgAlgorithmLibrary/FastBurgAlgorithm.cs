using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    ///     Fast implimentation of Burg algorithm for real signals.
    ///     For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    ///     FastBurgAlgorithm uses internal variables of type double
    /// </summary>
    public class FastBurgAlgorithm<T> where T : IConvertible
    {
        private readonly T[] _xInputSignal;

        /// <summary>
        ///     Position in x_inputSignal that we need prediction for.
        /// </summary>
        private int _absolutePosition;

        private T[] _aPredictionCoefs;
        private T[] _c;

        /// <summary>
        ///     Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        private T[] _deltaRAndAProduct;

        private T[] _g;

        // Naming: 
        // - first letter is the same as in the Koen Vos paper
        // - the part after underscore is my description
        private int _iIterationCounter;
        private T[] _kReflectionCoefs;
        private int _mCoefficientsNumber;
        private int _nHistoryLengthSamples;
        private T[] _r;
        private IOperations<T> _operations;

        public FastBurgAlgorithm(double[] inputSignal)
        {
            _xInputSignal = new T[inputSignal.Length];
            for (var index = 0; index < inputSignal.Length; index++)
                _xInputSignal[index] = 
                    (T)Convert.ChangeType(inputSignal[index], typeof(T));
            _operations = Operations<T>.Default;
        }

        /// <summary>
        ///     Calculates prediction coefficients for one sample using CPU
        /// </summary>
        /// <param name="position">
        ///     Position in inputSignal that we need
        ///     prediction for. Must be greater than historyLengthSamples
        /// </param>
        /// <param name="coefficientsNumber">
        ///     Number of prediction coefficients
        ///     that will be calculated. Greater number gives more accurate
        ///     prediction but takes more time to calculate
        /// </param>
        /// <param name="historyLengthSamples">
        ///     Number of samples that will
        ///     be used to calculate prediction coefficients
        /// </param>
        public void Train(
            int position,
            int coefficientsNumber,
            int historyLengthSamples)
        {
            _absolutePosition = position;
            _mCoefficientsNumber = coefficientsNumber;
            _nHistoryLengthSamples = historyLengthSamples;

            CreateInternalVariables();

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
        ///     Creates internal variables with desirable length
        /// </summary>
        private void CreateInternalVariables()
        {
            _aPredictionCoefs = new T[_mCoefficientsNumber + 1];
            _g = new T[_mCoefficientsNumber + 2];
            _r = new T[_mCoefficientsNumber + 1];
            _c = new T[_mCoefficientsNumber + 1];
            _kReflectionCoefs = new T[_mCoefficientsNumber];
            _deltaRAndAProduct = new T[_mCoefficientsNumber + 1];
        }

        /// <summary>
        ///     Returns forward prediction based on prediction coefficients that were
        ///     previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double GetForwardPrediction()
        {
            double prediction = 0;
            for (var index = 1; index <= _aPredictionCoefs.Length - 1; index++)
            {
                var multResult = _operations.Multiply(
                    _aPredictionCoefs[index],
                    _xInputSignal[_absolutePosition - index]);

                prediction -= (double) Convert.ChangeType(
                    multResult,
                    typeof(double));
            }

            return prediction;
        }

        /// <summary>
        ///     Returns backward prediction based on prediction coefficients that were
        ///     previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double GetBackwardPrediction()
        {
            double prediction = 0;
            for (var index = 1; index <= _aPredictionCoefs.Length - 1; index++)
            {
                var multResult = _operations.Multiply(
                    _aPredictionCoefs[index],
                    _xInputSignal[_absolutePosition -
                                  _nHistoryLengthSamples - 1 +
                                  index]);
                prediction -= (double)Convert.ChangeType(
                    multResult,
                    typeof(double));
            }

            return prediction;
        }

        /// <summary>
        ///     Returns prediction coefficients that were
        ///     previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public T[] GetPredictionCoefs()
        {
            var predictionCoefs = (T[]) _aPredictionCoefs.Clone();

            return predictionCoefs;
        }

        /// <summary>
        ///     Returns prediction coefficients that were
        ///     previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public T[] GetReflectionCoefs()
        {
            var reflectionCoefs = (T[]) _kReflectionCoefs.Clone();

            return reflectionCoefs;
        }

        /// <summary>
        ///     Updates vector g. For details see step 7 of algorithm on page 3 of
        ///     A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateG()
        {
            var oldG = (T[]) _g.Clone();

            // g.Length is i_iterationCounter + 1
            for (var index = 0; index <= _iIterationCounter; index++)
                _g[index] =
                    _operations.Add(
                        _operations.Add(
                            oldG[index],
                            _operations.Multiply(
                                _kReflectionCoefs[_iIterationCounter - 1],
                                oldG[JinversOrder(index, _iIterationCounter)])),
                        _deltaRAndAProduct[index]);

            for (var index = 0; index <= _iIterationCounter; index++)
                _g[_iIterationCounter + 1] = 
                    _operations.Add(
                        _g[_iIterationCounter + 1],
                        _operations.Multiply(
                            _r[index],
                            _aPredictionCoefs[index]));
        }

        /// <summary>
        ///     Calculates vector deltaRAndAProduct. For details see step 6 of algorithm on page 3 of
        ///     A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void ComputeDeltaRMultByA()
        {
            for (var indexRow = 0; indexRow <= _iIterationCounter; indexRow++)
            {
                var innerProduct1 = (T) Convert.ChangeType(0, typeof(T));
                var innerProduct2 = (T) Convert.ChangeType(0, typeof(T));
                for (var indexColumn = 0;
                    indexColumn <= _iIterationCounter;
                    indexColumn++)
                {
                    innerProduct1 =
                        _operations.Add(
                            innerProduct1,
                            _operations.Multiply(
                                _xInputSignal[_absolutePosition - _nHistoryLengthSamples +
                                              _iIterationCounter - indexColumn],
                                _aPredictionCoefs[indexColumn]));
                    innerProduct2 =
                        _operations.Add(
                            innerProduct2,
                            _operations.Multiply(
                                _xInputSignal[_absolutePosition - 1 -
                                      _iIterationCounter + indexColumn],
                                _aPredictionCoefs[indexColumn]));
                }

                _deltaRAndAProduct[indexRow] =
                    _operations.Subtract(
                        _operations.Subtract(
                            (T) Convert.ChangeType(0, typeof(T)),
                            _operations.Multiply(
                                _xInputSignal[_absolutePosition - _nHistoryLengthSamples +
                                   _iIterationCounter - indexRow],
                                innerProduct1)),
                        _operations.Multiply(
                            _xInputSignal[_absolutePosition - 1 -
                                  _iIterationCounter + indexRow],
                            innerProduct2));
            }
        }

        /// <summary>
        ///     Updates vector r. For details see step 5 of algorithm on page 3 of
        ///     A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdateR()
        {
            var oldR = (T[]) _r.Clone();

            for (var index = 0; index <= _iIterationCounter - 1; index++)
                _r[index + 1] = _operations.Subtract(
                                    _operations.Subtract(
                                        oldR[index],
                                        _operations.Multiply(
                                            _xInputSignal[_absolutePosition - _nHistoryLengthSamples + index],
                                            _xInputSignal[_absolutePosition - _nHistoryLengthSamples + _iIterationCounter])),
                                    _operations.Multiply(
                                        _xInputSignal[_absolutePosition - 1 - index],
                                        _xInputSignal[_absolutePosition - 1 - _iIterationCounter]));

            _r[0] = _operations.Multiply(
                        (T) Convert.ChangeType(2, typeof(T)),
                        _c[_iIterationCounter + 1]);
        }

        /// <summary>
        ///     Updates vector of prediction coefficients. For details see step 2 of
        ///     algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void UpdatePredictionCoefs()
        {
            var oldAPredictionCoefs = (T[]) _aPredictionCoefs.Clone();

            for (var index = 0; index <= _iIterationCounter + 1; index++)
                _aPredictionCoefs[index] = _operations.Add(
                                               oldAPredictionCoefs[index],
                                               _operations.Multiply(
                                                   _kReflectionCoefs[_iIterationCounter],
                                                   oldAPredictionCoefs[JinversOrder(index, _iIterationCounter + 1)]));
        }

        /// <summary>
        ///     Computes vector of reflection coefficients. For details see step 1
        ///     of algorithm on page 3 of A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void ComputeReflectionCoef()
        {
            var nominator = (T) Convert.ChangeType(0, typeof(T));
            var denominator = _operations.AlmostZero();

            for (var index = 0; index <= _iIterationCounter + 1; index++)
            {
                nominator = _operations.Add(
                                nominator,
                                _operations.Multiply(
                                    _aPredictionCoefs[index],
                                     _g[JinversOrder(index, _iIterationCounter + 1)]));
                denominator = _operations.Add(
                                    denominator,
                                    _operations.Multiply(
                                        _aPredictionCoefs[index],
                                        _g[index]));
            }

            _kReflectionCoefs[_iIterationCounter] = _operations.Divide(
                            _operations.Subtract(
                                (T) Convert.ChangeType(0, typeof(T)), 
                                nominator),
                            denominator);
        }

        /// <summary>
        ///     Inverts index to flip a vector insted of multiplication with J matrix.
        ///     For details see (12) on page 2 of
        ///     A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        /// <param name="index">from 0 to max</param>
        /// <param name="max">positive number</param>
        /// <returns></returns>
        private static int JinversOrder(int index, int max)
        {
            return max - index;
        }

        /// <summary>
        ///     Initializes i_iterationCounter and vectors. For details see step 0 of
        ///     algorithm on page 3 of
        ///     A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void Initialization()
        {
            FindAutocorrelations();

            _iIterationCounter = 0;
            _aPredictionCoefs[0] = (T) Convert.ChangeType(1, typeof(T));
            _g[0] = _operations.Subtract(
                        _operations.Subtract(
                            _operations.Multiply(
                                (T) Convert.ChangeType(2, typeof(T)),
                                _c[0]),
                            _operations.Multiply(
                                _operations.Abs(_xInputSignal[_absolutePosition - _nHistoryLengthSamples]),
                                _operations.Abs(_xInputSignal[_absolutePosition - _nHistoryLengthSamples]))),
                        _operations.Multiply(
                            _operations.Abs(_xInputSignal[_absolutePosition - 1]),
                            _operations.Abs(_xInputSignal[_absolutePosition - 1])));
            _g[1] = _operations.Multiply(
                        (T)Convert.ChangeType(2, typeof(T)),
                        _c[1]);
            // the paper says r[1], error in paper?
            _r[0] = _operations.Multiply(
                        (T)Convert.ChangeType(2, typeof(T)),
                        _c[1]);
        }

        /// <summary>
        ///     Calculates autocorrelations. For details see step 0 of
        ///     algorithm on page 3 of
        ///     A Fast Implementation of Burg’s Method by Koen Vos
        /// </summary>
        private void FindAutocorrelations()
        {
            for (var j = 0; j <= _mCoefficientsNumber; j++)
            {
                _c[j] = (T)Convert.ChangeType(0, typeof(T));
                for (var index = _absolutePosition - _nHistoryLengthSamples;
                    index <= _absolutePosition - 1 - j;
                    index++)
                    _c[j] = _operations.Add(
                                _c[j],
                                _operations.Multiply(
                                    _xInputSignal[index],
                                    _xInputSignal[index + j]));
            }
        }
    }
}