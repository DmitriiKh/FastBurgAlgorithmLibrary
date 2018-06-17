using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Fast implimentation of Burg algorithm for real signals.
    /// For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    /// FastBurgAlgorithm64 uses internal variables of type double
    /// </summary>
    public class FastBurgAlgorithm64 : FastBurgAlgorithm
    {
        public FastBurgAlgorithm64(double[] inputSignal) :
            base(inputSignal, typeof(double))
        {

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
        /// Returns reflection coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public double[] GetReflectionCoefs()
        {
            double[] reflectionCoefs = (double[])k_reflectionCoefs.Clone();

            return reflectionCoefs;
        }
    }
}
