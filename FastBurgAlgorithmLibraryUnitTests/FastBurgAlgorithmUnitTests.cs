using System;
using NUnit.Framework;
using FastBurgAlgorithmLibrary;

namespace FastBurgAlgorithmLibraryUnitTests
{
    [TestFixture]
    public class FastBurgAlgorithmUnitTests
    {
        [Test]
        [TestCase(64, 0.00001)]
        [TestCase(128, 0.0000001)]
        public void GetForwardPrediction_SinInput_ReturnsCorrectPrediction(int precision, double accuracy)
        {
            const int coefNumber = 4;
            const int historyLength = 512;
            const int numberOfSamplesToCheck = 10;

            var inputAudio = 
                new double[historyLength + numberOfSamplesToCheck]; 

            for (var i = 0; i < inputAudio.Length; i++)
            {
                inputAudio[i] = Math.Sin( 
                    2 * Math.PI * i / (historyLength / 5.2));
            }

            var fba = new FastBurgAlgorithm(inputAudio);

            for (var index = historyLength;
                index < historyLength + numberOfSamplesToCheck;
                index++)
            {
                fba.Train(precision, index, coefNumber, historyLength);
                var forwardPrediction = fba.GetForwardPrediction();

                Assert.AreEqual(
                    inputAudio[index],
                    forwardPrediction,
                    accuracy);
            }
        }
    }
}
