using NUnit.Framework;
using FastBurgAlgorithmLibrary;

namespace FastBurgAlgorithmLibraryUnitTests
{
    [TestFixture]
    public class FastBurgAlgorithm64UnitTests
    {
        [Test]
        public void GetForwardPrediction_SinInput_ReturnsCorrectPrediction()
        {
            const int coefNumber = 4;
            const int historyLength = 512;
            const int numberOfSamplesToCheck = 10;

            double[] inputAudio = 
                new double[historyLength + numberOfSamplesToCheck]; 

            for (int i = 0; i < inputAudio.Length; i++)
            {
                inputAudio[i] = System.Math.Sin( 
                    2 * System.Math.PI * i / (historyLength / 5.2));
            }

            FastBurgAlgorithm64 fba = new FastBurgAlgorithm64(inputAudio);

            for (int index = historyLength;
                index < historyLength + numberOfSamplesToCheck;
                index++)
            {
                fba.Train(index, coefNumber, historyLength);
                var forwardPrediction = fba.GetForwardPrediction();

                Assert.AreEqual(
                    inputAudio[index],
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

            const double accuracy = 0.01;

            double[] inputAudio = 
                new double[historyLength + numberOfSamplesToCheck]; 

            for (int i = 0; i < inputAudio.Length; i++)
            {
                inputAudio[i] = System.Math.Sin( 
                    2 * System.Math.PI * i / (historyLength / 5.2));
            }

            FastBurgAlgorithm64 fba = new FastBurgAlgorithm64(inputAudio);

            fba.Train(historyLength, coefNumber, historyLength);

            double[] predictionCoefs = fba.GetPredictionCoefs();

            Assert.AreEqual(
                1,
                predictionCoefs[0],
                accuracy);

            Assert.AreEqual(
                -3.991510267867756,
                predictionCoefs[1],
                accuracy);

            Assert.AreEqual(
                5.983035128379795,
                predictionCoefs[2],
                accuracy);

            Assert.AreEqual(
                -3.991503459864878,
                predictionCoefs[3],
                accuracy);

            Assert.AreEqual(
                0.9999965889050035,
                predictionCoefs[4],
                accuracy);

            double[] reflectionCoefs = fba.GetReflectionCoefs();

            Assert.AreEqual(
                -0.9979213453536945,
                reflectionCoefs[0],
                accuracy);

            Assert.AreEqual(
                0.9999990984096440,
                reflectionCoefs[1],
                accuracy);

            Assert.AreEqual(
                -0.9978363901155060,
                reflectionCoefs[2],
                accuracy);

            Assert.AreEqual(
                0.9999965889050035,
                reflectionCoefs[3],
                accuracy);
        }

        [Test]
        public void GetForwardPrediction_ZeroInput_ReturnsCorrectPrediction()
        {
            const int coefNumber = 4;
            const int historyLength = 512;
            const int numberOfSamplesToCheck = 10;

            double[] inputAudio =
                new double[historyLength + numberOfSamplesToCheck];

            for (int i = 0; i < inputAudio.Length; i++)
            {
                inputAudio[i] = 0;
            }

            FastBurgAlgorithm64 fba = new FastBurgAlgorithm64(inputAudio);

            for (int index = historyLength;
                index < historyLength + numberOfSamplesToCheck;
                index++)
            {
                fba.Train(index, coefNumber, historyLength);
                var forwardPrediction = fba.GetForwardPrediction();

                Assert.AreEqual(
                    inputAudio[index],
                    forwardPrediction,
                    0.0000001);
            }
        }
    }
}
