using NUnit.Framework;
using FastBurgAlgorithmLibrary;

namespace FastBurgAlgorithmLibraryUnitTests
{
    [TestFixture]
    public class FastBurgAlgorithmUnitTests
    {
        [Test]
        public void GetForwardPrediction_SinInput_ReturnsCorrectPrediction()
        {
            const int coefNumber = 4;
            const int historyLength = 512;
            const int numberOfSamplesToCheck = 10;

            double[] input_audio = // float
                new double[historyLength + numberOfSamplesToCheck]; // float

            for (int i = 0; i < input_audio.Length; i++)
            {
                input_audio[i] = System.Math.Sin( // float
                    2 * System.Math.PI * i / (historyLength / 5.2));
            }

            FastBurgAlgorithm fba = new FastBurgAlgorithm(input_audio);

            for (int index = historyLength;
                index < historyLength + numberOfSamplesToCheck;
                index++)
            {
                fba.Train(index, coefNumber, historyLength);
                var forwardPrediction = fba.GetForwardPrediction();

                Assert.AreEqual(
                    input_audio[index],
                    forwardPrediction,
                    0.0000001);
            }
        }

        [Test]
        public void Train_SinInput_ReturnsCorrectPredictionAndReflectionCoefficients()
        {
            /* Coefficients for comparision are taken from GNU Octave arburg() function
             * t = [0:2000]
             * x = sin( 2 * pi() * t / (512 / 5.2))
             * output_precision(16)
             * [a, v, k] = arburg(x(1:512),4)
             */
            const int coefNumber = 4;
            const int historyLength = 512;
            const int numberOfSamplesToCheck = 1;

            double[] input_audio = 
                new double[historyLength + numberOfSamplesToCheck]; 

            for (int i = 0; i < input_audio.Length; i++)
            {
                input_audio[i] = System.Math.Sin( 
                    2 * System.Math.PI * i / (historyLength / 5.2));
            }

            FastBurgAlgorithm fba = new FastBurgAlgorithm(input_audio);

            fba.Train(historyLength, coefNumber, historyLength);

            decimal[] predictionCoefs = fba.GetPredictionCoefs();

            Assert.AreEqual(
                1,
                (double)predictionCoefs[0],
                0.00001);

            Assert.AreEqual(
                -3.991510267867756,
                (double)predictionCoefs[1],
                0.00001);

            Assert.AreEqual(
                5.983035128379795,
                (double)predictionCoefs[2],
                0.00001);

            Assert.AreEqual(
                -3.991503459864878,
                (double)predictionCoefs[3],
                0.00001);

            Assert.AreEqual(
                0.9999965889050035,
                (double)predictionCoefs[4],
                0.00001);

            decimal[] reflectionCoefs = fba.GetReflectionCoefs();

            Assert.AreEqual(
                -0.9979213453536945,
                (double)reflectionCoefs[0],
                0.000001);

            Assert.AreEqual(
                0.9999990984096440,
                (double)reflectionCoefs[1],
                0.00001);

            Assert.AreEqual(
                -0.9978363901155060,
                (double)reflectionCoefs[2],
                0.00001);

            Assert.AreEqual(
                0.9999965889050035,
                (double)reflectionCoefs[3],
                0.00001);
        }
    }
}
