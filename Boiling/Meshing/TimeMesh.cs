namespace Boiling.Meshing;

public class TimeMesh
{
    private readonly double[] _times;
    public double this[int timeMoment] => _times[timeMoment];

    public int TimesCount => _times.Length;

    public TimeMesh(double start, double end, int splits)
    {
        double step = (end - start) / splits;

        _times = new double[splits + 1];

        for (int i = 0; i < splits + 1; i++)
        {
            _times[i] = start + i * step;
        }
    }
}