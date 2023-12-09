using System.Collections;

namespace Boiling.MathHelper;

public class Vector(int length) : IEnumerable<double>
{
    private readonly double[] _storage = new double[length];
    public int Length { get; } = length;

    public double this[int idx]
    {
        get => _storage[idx];
        set => _storage[idx] = value;
    }

    public void FromCollection(IEnumerable<double> collection)
    {
        int idx = 0;
        foreach (var value in collection)
        {
            _storage[idx++] = value;
        }
    }

    public static double Dot(Vector a, Vector b)
        => a.Select((t, i) => t * b[i]).Sum();

    public static double Dot(double[] a, double[] b)
        => a.Select((t, i) => t * b[i]).Sum();

    public static void Copy(Vector source, Vector? destination)
    {
        destination ??= new Vector(source.Length);
        
        for (int i = 0; i < source.Length; i++)
        {
            destination[i] = source[i];
        }
    }

    public void Fill(double value = 0.0)
        => Array.Fill(_storage, value);

    public double Norm()
    {
        double result = 0.0;

        for (int i = 0; i < Length; i++)
        {
            result += _storage[i] * _storage[i];
        }

        return Math.Sqrt(Convert.ToDouble(result));
    }

    public double SqrNorm()
    {
        double result = 0.0;

        for (int i = 0; i < Length; i++)
        {
            result += _storage[i] * _storage[i];
        }

        return result;
    }
    
    public IEnumerator<double> GetEnumerator()
        => ((IEnumerable<double>)_storage).GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}