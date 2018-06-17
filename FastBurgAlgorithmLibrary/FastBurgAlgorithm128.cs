using System;

namespace FastBurgAlgorithmLibrary
{
    /// <summary>
    /// Fast implimentation of Burg algorithm for real signals.
    /// For details see paper A Fast Implementation of Burg’s Method by Koen Vos.
    /// FastBurgAlgorithm128 uses internal variables of type decimal
    /// </summary>
    public class FastBurgAlgorithm128 : FastBurgAlgorithm
    {
        public FastBurgAlgorithm128(double[] inputSignal): 
            base(inputSignal, typeof(decimal))
        {
            
        }

        /// <summary>
        /// Returns prediction coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public decimal[] GetPredictionCoefs()
        {
            decimal[] predictionCoefs = (decimal[])a_predictionCoefs.Clone();

            return predictionCoefs;
        }

        /// <summary>
        /// Returns reflection coefficients that were
        /// previously calculated with Train() method
        /// </summary>
        /// <returns></returns>
        public decimal[] GetReflectionCoefs()
        {
            decimal[] reflectionCoefs = (decimal[])k_reflectionCoefs.Clone();

            return reflectionCoefs;
        }
    }
}
