using System;

namespace FastBurgAlgorithmLibrary
{
    internal abstract class Data
    {
        internal dynamic g;
        internal dynamic old_g;
        internal dynamic r;
        internal dynamic old_r;
        internal dynamic c;
        internal dynamic k_reflectionCoefs;
        internal dynamic a_predictionCoefs;
        internal dynamic old_a_predictionCoefs;
        internal Type InternalVariableType;

        /// <summary>
        /// Product of deltaR matrix and a_predictionCoefs
        /// </summary>
        internal dynamic deltaRAndAProduct;

        internal abstract dynamic DoubleToInternal(double value);

        internal double InternalToDouble(ValueType value)
        {
            return (double) value;
        }

        internal abstract ValueType Abs(ValueType value);

        internal abstract ValueType Multiply(ValueType operand1, ValueType operand2);
    }
}
