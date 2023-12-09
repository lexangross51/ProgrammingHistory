namespace Boiling.MathHelper;

public abstract class IterativeSolver
{
    protected SparseMatrix Matrix = default!;
    protected double[] RightPart = default!;
    protected readonly int MaxIterations;
    protected readonly double Eps;
    protected double[]? YVector;
    
    public int IterationsCount { get; protected set; }
    public double[]? Solution { get; protected set; }

    protected IterativeSolver(int maxIterations, double eps)
        => (MaxIterations, Eps) = (maxIterations, eps);

    public void SetSystem(SparseMatrix matrix, double[] rightPart)
        => (Matrix, RightPart) = (matrix, rightPart);

    public abstract void Compute();

    protected void CholeskyDecomposition(double[] diNew, double[] ggNew)
    {
        double sumLower = 0.0;
        double sumDi = 0.0;

        for (int i = 0; i < Matrix.Size; i++)
        {
            int i0 = Matrix.Ig[i];
            int i1 = Matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
            {
                int j = Matrix.Jg[k];
                int j0 = Matrix.Ig[j];
                int j1 = Matrix.Ig[j + 1];
                int ik = i0;
                int kj = j0;

                while (ik < k && kj < j1)
                {
                    if (Matrix.Jg[ik] == Matrix.Jg[kj])
                    {
                        sumLower += ggNew[ik] * ggNew[kj];
                        ik++;
                        kj++;
                    }
                    else
                    {
                        if (Matrix.Jg[ik] > Matrix.Jg[kj])
                            kj++;
                        else
                            ik++;
                    }
                }

                ggNew[k] = (ggNew[k] - sumLower) / diNew[j];
                sumDi += ggNew[k] * ggNew[k];
                sumLower = 0.0;
            }

            diNew[i] = Math.Sqrt(diNew[i] - sumDi);
            sumDi = 0.0;
        }
    }

    protected void MoveForCholesky(double[] diNew, double[] ggNew, double[] vector, double[] result)
    {
        YVector ??= new double[RightPart.Length];
        Array.Copy(vector, YVector, vector.Length);
        double sum = 0.0;

        for (int i = 0; i < Matrix.Size; i++)
        {
            int i0 = Matrix.Ig[i];
            int i1 = Matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
                sum += ggNew[k] * YVector[Matrix.Jg[k]];

            YVector[i] = (YVector[i] - sum) / diNew[i];
            sum = 0.0;
        }

        Array.Copy(YVector, result, YVector.Length);

        for (int i = Matrix.Size - 1; i >= 0; i--)
        {
            int i0 = Matrix.Ig[i];
            int i1 = Matrix.Ig[i + 1];
            result[i] = YVector[i] / diNew[i];

            for (int k = i0; k < i1; k++)
                YVector[Matrix.Jg[k]] -= ggNew[k] * result[i];
        }
    }
}

public class LOS(int maxIterations, double eps) 
    : IterativeSolver(maxIterations, eps)
{
    public override void Compute()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(Matrix, $"{nameof(Matrix)} cannot be null, set the matrix");
            ArgumentNullException.ThrowIfNull(RightPart, $"{nameof(RightPart)} cannot be null, set the vector");

            Solution = new double[RightPart.Length];

            double[] z = new double[RightPart.Length];
            double[] r = new double[RightPart.Length];
            double[] p = new double[RightPart.Length];
            double[] product = new double[RightPart.Length];
            SparseMatrix.Dot(Matrix, Solution, product);

            for (int i = 0; i < product.Length; i++)
            {
                r[i] = RightPart[i] - product[i];
            }

            Array.Copy(r, z, r.Length);
            SparseMatrix.Dot(Matrix, z, p);

            var squareNorm = Vector.Dot(r, r);

            for (IterationsCount = 0; IterationsCount < MaxIterations && squareNorm > Eps; IterationsCount++)
            {
                var alpha = Vector.Dot(p, r) / Vector.Dot(p, p);

                for (int i = 0; i < Solution.Length; i++)
                {
                    Solution[i] += alpha * z[i];
                }

                squareNorm = Vector.Dot(r, r) - (alpha * alpha * Vector.Dot(p, p));

                for (int i = 0; i < Solution.Length; i++)
                {
                    r[i] -= alpha * p[i];
                }

                SparseMatrix.Dot(Matrix, r, product);

                var beta = -Vector.Dot(p, product) / Vector.Dot(p, p);

                for (int i = 0; i < Solution.Length; i++)
                {
                    z[i] = r[i] + beta * z[i];
                    p[i] = product[i] + beta * p[i];
                }
            }
        }
        catch (ArgumentNullException ex)
        {
            Console.WriteLine($"We had problem: {ex.Message}");
            throw;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"We had problem: {ex.Message}");
        }
    }
}

public class CGMCholesky(int maxIterations, double eps) 
    : IterativeSolver(maxIterations, eps)
{
    public override void Compute()
    {
        ArgumentNullException.ThrowIfNull(Matrix, $"{nameof(Matrix)} cannot be null, set the matrix");
        ArgumentNullException.ThrowIfNull(RightPart, $"{nameof(Vector)} cannot be null, set the vector");

        double rightPartNorm = RightPart.Norm();
        Solution = new double[RightPart.Length];
        Array.Fill(Solution, 1.0);
        var r = new double[RightPart.Length];
        var z = new double[RightPart.Length];
        var result = new double[RightPart.Length];
        var product = new double[RightPart.Length];
        double[] diNew = new double[Matrix.Size];
        double[] ggNew = new double[Matrix.Ggl.Length];
        Array.Copy(Matrix.Di, diNew, Matrix.Size);
        Array.Copy(Matrix.Ggl, ggNew, Matrix.Ggl.Length);
        
        CholeskyDecomposition(diNew, ggNew);
        
        SparseMatrix.Dot(Matrix, Solution, product);

        for (int i = 0; i < r.Length; i++)
        {
            r[i] = RightPart[i] - product[i];
        }
        
        MoveForCholesky(diNew, ggNew, r, z);

        for (IterationsCount = 0; IterationsCount < MaxIterations && r.Norm() / rightPartNorm >= Eps; IterationsCount++)
        {
            MoveForCholesky(diNew, ggNew, r, result);
            
            SparseMatrix.Dot(Matrix, z, product);
            double mrDotR = Vector.Dot(result, r);
            double alpha = mrDotR / Vector.Dot(product, z);

            for (int i = 0; i < Solution.Length; i++)
            {
                Solution[i] += alpha * z[i];
                r[i] -= alpha * product[i];
            }
            
            MoveForCholesky(diNew, ggNew, r, result);

            double beta = Vector.Dot(result, r) / mrDotR;

            for (int i = 0; i < z.Length; i++)
            {
                z[i] = result[i] + beta * z[i];
            }
        }
    }
}