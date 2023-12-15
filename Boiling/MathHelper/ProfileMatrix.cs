namespace Boiling.MathHelper;

public class ProfileMatrix
{
    public int[] Ig { get; init; } = null!;
    public double[] Di { get; init; } = null!;
    public double[] Ggl { get; init; } = null!;
    public double[] Ggu { get; init; } = null!;
    public int Size => Di.Length;
    
    public static void Dot(ProfileMatrix matrix, double[] vector, double[]? product)
    {
        if (matrix.Size != vector.Length)
        {
            throw new Exception("Size of matrix not equal to size of vector");
        }

        product ??= new double[vector.Length];
        Array.Fill(product, 0.0);
        
        for (int i = 0; i < product.Length; i++)
        {
            product[i] = matrix.Di[i] * vector[i];

            int l = matrix.Ig[i + 1] - matrix.Ig[i];
            int k = i - 1;

            for (int j = 0; j < l; j++)
            {
                int index = matrix.Ig[i] + j - 1;

                product[i] += matrix.Ggl[index] * vector[k];
                product[k] += matrix.Ggu[index] * vector[i];
            }
        }
    }
}