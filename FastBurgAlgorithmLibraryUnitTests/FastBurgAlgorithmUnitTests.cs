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

            float[] input_audio =
                new float[historyLength + numberOfSamplesToCheck];

            for (int i = 0; i < input_audio.Length; i++)
            {
                input_audio[i] = (float)System.Math.Sin(
                    2 * System.Math.PI * i / (historyLength / 5.2));
            }

            for (int index = historyLength;
                index < historyLength + numberOfSamplesToCheck;
                index++)
            {
                FastBurgAlgorithm fbp = new FastBurgAlgorithm(input_audio);
                fbp.Train(index, coefNumber, historyLength);
                var forwardPrediction = fbp.GetForwardPrediction();

                Assert.AreEqual(
                    input_audio[index],
                    forwardPrediction,
                    0.000001);
            }
        }
    }
}
