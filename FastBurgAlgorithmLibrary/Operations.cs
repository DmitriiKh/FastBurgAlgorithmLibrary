using System;

namespace FastBurgAlgorithmLibrary
{
    public interface IOperations<T>
    {
        T Add(T a, T b);
        T Subtract(T a, T b);
        T Multiply(T a, T b);
        T Divide(T a, T b);
        T Abs(T a);
        T AlmostZero();
    }

    public static class Operations<T>
    {
        public static IOperations<T> Default => Create();

        static IOperations<T> Create()
        {
            var type = typeof(T);
            switch (Type.GetTypeCode(type))
            {
                case TypeCode.Double:
                    return (IOperations<T>)new DoubleOperations();
                case TypeCode.Decimal:
                    return (IOperations<T>)new DecimalOperations();
                default:
                    var message = $"Operations for type {type.Name} is not supported.";
                    throw new NotSupportedException(message);
            }
        }

        class DoubleOperations : IOperations<double>
        {
            public double Add(double a, double b) => a + b;
            public double Subtract(double a, double b) => a - b;
            public double Multiply(double a, double b) => a * b;
            public double Divide(double a, double b) => a / b;
            public double Abs(double a) => Math.Abs(a);
            public double AlmostZero() => double.Epsilon;
        }

        class DecimalOperations : IOperations<decimal>
        {
            public decimal Add(decimal a, decimal b) => a + b;
            public decimal Subtract(decimal a, decimal b) => a - b;
            public decimal Multiply(decimal a, decimal b) => a * b;
            public decimal Divide(decimal a, decimal b) => a / b;
            public decimal Abs(decimal a) => Math.Abs(a);
            public decimal AlmostZero() => 0.000000000000000000000000001M;
        }
    }
}
