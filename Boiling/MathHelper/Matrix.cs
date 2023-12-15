namespace Boiling.MathHelper;

public class Matrix
{
    private readonly double[][] _storage;
    public int Rows { get; }
    public int Columns { get; }
    
    public Matrix(int nRows, int nColumns)
    {
        Rows = nRows;
        Columns = nColumns;
        _storage = new double[Rows].Select(_ => new double[nColumns]).ToArray();
    }
    
    public double this[int i, int j]
    {
        get => _storage[i][j];
        set => _storage[i][j] = value;
    }

    public void Fill(double value)
    {
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Columns; j++)
            {
                _storage[i][j] = value;
            }
        }
    }
    
    public static void Dot(Matrix matrix, IEnumerable<double> vector, double[]? product)
    {
        var v = vector.ToArray();
        
        if (matrix.Columns != v.Length)
        {
            throw new Exception("Numbers of columns not equal to size of vector");
        }

        product ??= new double[v.Length];
        Array.Fill(product, 0.0);

        for (int i = 0; i < matrix.Rows; i++)
        {
            for (int j = 0; j < matrix.Columns; j++)
            {
                product[i] += matrix[i, j] * v[j];
            }
        }
    }
}