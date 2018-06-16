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
                    0.0001);
            }
        }
    }
}
