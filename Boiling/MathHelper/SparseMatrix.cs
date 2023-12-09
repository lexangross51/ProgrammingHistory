namespace Boiling.MathHelper;

public class SparseMatrix
{
    public int[] Ig { get; init; }
    public int[] Jg { get; init; }
    public double[] Di { get; }
    public double[] Ggl { get; }
    public double[] Ggu { get; }
    public int Size => Di.Length;
    
    public SparseMatrix(int[] ig, int[] jg)
    {
        Ig = ig;
        Jg = jg;
        Di = new double[Ig.Length - 1];
        Ggl = new double[Jg.Length];
        Ggu = new double[Jg.Length];
    }
    
    public void Add(int i, int j, double value)
    {
        if (i == j)
        {
            Di[i] += value;
        }
        else if (i < j)
        {
            for (int idx = Ig[j]; idx < Ig[j + 1]; idx++)
            {
                if (Jg[idx] != i) continue;
                Ggu[idx] += value;
            }
        }
        else
        {
            for (int idx = Ig[i]; idx < Ig[i + 1]; idx++)
            {
                if (Jg[idx] != j) continue;
                Ggl[idx] += value;
            }
        }
    }
    
    public static void Dot(SparseMatrix matrix, double[] vector, double[]? product)
    {
        if (matrix.Size != vector.Length)
        {
            throw new Exception("Size of matrix not equal to size of vector");
        }

        product ??= new double[vector.Length];
        Array.Fill(product, 0.0);
        var ig = matrix.Ig;
        var jg = matrix.Jg;
        var di = matrix.Di;
        var ggl = matrix.Ggl;
        var ggu = matrix.Ggu;

        for (int i = 0; i < vector.Length; i++)
        {
            product[i] = di[i] * vector[i];

            for (int j = ig[i]; j < ig[i + 1]; j++)
            {
                product[i] += ggl[j] * vector[jg[j]];
                product[jg[j]] += ggu[j] * vector[i];
            }
        }
    }
    
    public void PrintDense(string path)
    {
        double[,] a = new double[Size, Size];

        for (int i = 0; i < Size; i++)
        {
            a[i, i] = Di[i];

            for (int j = Ig[i]; j < Ig[i + 1]; j++)
            {
                a[i, Jg[j]] = Ggl[j];
                a[Jg[j], i] = Ggu[j];
            }
        }

        using var sw = new StreamWriter(path);
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size; j++)
            {
                sw.Write(a[i, j].ToString("0.0000000") + "\t");
            }

            sw.WriteLine();
        }
    }

    public void Clear()
    {
        for (int i = 0; i < Size; i++)
        {
            Di[i] = 0.0;
        }

        for (int i = 0; i < Ggl.Length; i++)
        {
            Ggl[i] = 0.0;
            Ggu[i] = 0.0;
        }
    }
}