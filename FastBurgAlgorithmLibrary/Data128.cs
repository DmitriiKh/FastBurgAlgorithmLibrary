using System;

namespace FastBurgAlgorithmLibrary
{
    internal class Data128 : Data
    {
        public Data128(int mCoefficientsNumber)
        {
            var length = mCoefficientsNumber + 1;

            InternalVariableType = typeof(decimal);

            a_predictionCoefs = new decimal[length + 1];
            old_a_predictionCoefs = new decimal[length + 1];
            g = new decimal[length + 1];
            old_g = new decimal[length + 1];
            old_g = new decimal[length + 1];
            r = new decimal[length];
            old_r = new decimal[length];
            c = new decimal[length];
            k_reflectionCoefs = new decimal[length];
            deltaRAndAProduct = new decimal[length];
        }
        
        internal override ValueType Abs(ValueType value)
        {
            return Math.Abs((decimal)value);
        }

        internal override ValueType Multiply(ValueType operand1, ValueType operand2)
        {
            return (decimal) operand1 * (decimal) operand2;
        }

        internal override dynamic DoubleToInternal(double value)
        {
            return (decimal)value;
        }
    }
}
