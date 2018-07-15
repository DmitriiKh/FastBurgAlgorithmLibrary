using System;

namespace FastBurgAlgorithmLibrary
{
    internal class Data64 : Data
    {
        public Data64(int mCoefficientsNumber)
        {
            var length = mCoefficientsNumber + 1;

            InternalVariableType = typeof(double); 
            
            a_predictionCoefs = new double[length + 1];
            old_a_predictionCoefs = new double[length + 1];
            g = new double[length + 1];
            old_g = new double[length + 1];
            r = new double[length];
            old_r = new double[length];
            c = new double[length];
            k_reflectionCoefs = new double[length];
            deltaRAndAProduct = new double[length];
        }

        internal override ValueType Abs(ValueType value)
        {
            return Math.Abs((double) value);
        }

        internal override ValueType Multiply(ValueType operand1, ValueType operand2)
        {
            return (double) operand1 * (double) operand2;
        }

        internal override dynamic DoubleToInternal(double value)
        {
            return value;
        }
    }
}
