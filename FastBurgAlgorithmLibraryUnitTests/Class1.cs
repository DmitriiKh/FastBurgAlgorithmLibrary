using NUnit.Framework;

namespace FastBurgAlgorithmLibraryUnitTests
{
    [TestFixture]
    public class Class1
    {
        [Test]
        public void PassingTest()
        {
            Assert.AreEqual(4, Add(2, 2));
        }

        int Add(int x, int y)
        {
            return x + y;
        }
    }
}
