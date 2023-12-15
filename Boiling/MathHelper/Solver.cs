// namespace Boiling.MathHelper;
//
// public abstract class IterativeSolver
// {
//     protected SparseMatrix Matrix = default!;
//     protected double[] RightPart = default!;
//     protected readonly int MaxIterations;
//     protected readonly double Eps;
//     protected double[]? YVector;
//     
//     public int IterationsCount { get; protected set; }
//     public double[]? Solution { get; protected set; }
//
//     protected IterativeSolver(int maxIterations, double eps)
//         => (MaxIterations, Eps) = (maxIterations, eps);
//
//     public void SetSystem(SparseMatrix matrix, double[] rightPart)
//         => (Matrix, RightPart) = (matrix, rightPart);
//
//     public abstract void Compute();
//     
//     protected void LU(double[] gglnew, double[] ggunew, double[] dinew) {
//         double suml = 0.0;
//         double sumu = 0.0;
//         double sumdi = 0.0;
//
//         for (int i = 0; i < Matrix.Size; i++) {
//             int i0 = Matrix.Ig[i];
//             int i1 = Matrix.Ig[i + 1];
//
//             for (int k = i0; k < i1; k++) {
//                 int j = Matrix.Jg[k];
//                 int j0 = Matrix.Ig[j];
//                 int j1 = Matrix.Ig[j + 1];
//                 int ik = i0;
//                 int kj = j0;
//
//                 while (ik < k && kj < j1) {
//                     if (Matrix.Jg[ik] == Matrix.Jg[kj]) {
//                         suml += gglnew[ik] * ggunew[kj];
//                         sumu += ggunew[ik] * gglnew[kj];
//                         ik++;
//                         kj++;
//                     } else if (Matrix.Jg[ik] > Matrix.Jg[kj]) {
//                         kj++;
//                     } else {
//                         ik++;
//                     }
//                 }
//
//                 gglnew[k] -= suml;
//                 ggunew[k] = (ggunew[k] - sumu) / dinew[j];
//                 sumdi += gglnew[k] * ggunew[k];
//                 suml = 0.0;
//                 sumu = 0.0;
//             }
//
//             dinew[i] -= sumdi;
//             sumdi = 0.0;
//         }
//     }
//     
//     protected Vector Direct(Vector vector, double[] gglnew, double[] dinew) {
//         Vector y = new(vector.Length);
//         Vector.Copy(vector, y);
//
//         double sum = 0.0;
//
//         for (int i = 0; i < Matrix.Size; i++) {
//             int i0 = Matrix.Ig[i];
//             int i1 = Matrix.Ig[i + 1];
//
//             for (int k = i0; k < i1; k++)
//                 sum += gglnew[k] * y[Matrix.Jg[k]];
//
//             y[i] = (y[i] - sum) / dinew[i];
//             sum = 0.0;
//         }
//
//         return y;
//     }
//
//     protected Vector Reverse(Vector vector, double[] ggunew) {
//         Vector result = new(vector.Length);
//         Vector.Copy(vector, result);
//
//         for (int i = Matrix.Size - 1; i >= 0; i--) {
//             int i0 = Matrix.Ig[i];
//             int i1 = Matrix.Ig[i + 1];
//
//             for (int k = i0; k < i1; k++)
//                 result[Matrix.Jg[k]] -= ggunew[k] * result[i];
//         }
//
//         return result;
//     }
// }
//
// public class LOS(int maxIterations, double eps) 
//     : IterativeSolver(maxIterations, eps)
// {
//     public override void Compute()
//     {
//         try
//         {
//             ArgumentNullException.ThrowIfNull(Matrix, $"{nameof(Matrix)} cannot be null, set the matrix");
//             ArgumentNullException.ThrowIfNull(RightPart, $"{nameof(RightPart)} cannot be null, set the vector");
//
//             Solution = new double[RightPart.Length];
//
//             double[] z = new double[RightPart.Length];
//             double[] r = new double[RightPart.Length];
//             double[] p = new double[RightPart.Length];
//             double[] product = new double[RightPart.Length];
//             SparseMatrix.Dot(Matrix, Solution, product);
//
//             for (int i = 0; i < product.Length; i++)
//             {
//                 r[i] = RightPart[i] - product[i];
//             }
//
//             Array.Copy(r, z, r.Length);
//             SparseMatrix.Dot(Matrix, z, p);
//
//             var squareNorm = Vector.Dot(r, r);
//
//             for (IterationsCount = 0; IterationsCount < MaxIterations && squareNorm > Eps; IterationsCount++)
//             {
//                 var alpha = Vector.Dot(p, r) / Vector.Dot(p, p);
//
//                 for (int i = 0; i < Solution.Length; i++)
//                 {
//                     Solution[i] += alpha * z[i];
//                 }
//
//                 squareNorm = Vector.Dot(r, r) - (alpha * alpha * Vector.Dot(p, p));
//
//                 for (int i = 0; i < Solution.Length; i++)
//                 {
//                     r[i] -= alpha * p[i];
//                 }
//
//                 SparseMatrix.Dot(Matrix, r, product);
//
//                 var beta = -Vector.Dot(p, product) / Vector.Dot(p, p);
//
//                 for (int i = 0; i < Solution.Length; i++)
//                 {
//                     z[i] = r[i] + beta * z[i];
//                     p[i] = product[i] + beta * p[i];
//                 }
//             }
//         }
//         catch (ArgumentNullException ex)
//         {
//             Console.WriteLine($"We had problem: {ex.Message}");
//             throw;
//         }
//         catch (Exception ex)
//         {
//             Console.WriteLine($"We had problem: {ex.Message}");
//         }
//     }
// }
//
// // public class LOSLU(int maxIterations, double eps) : IterativeSolver(maxIterations, eps)
// // {
// //     public override void Compute()
// //     {
// //         try
// //         {
// //             ArgumentNullException.ThrowIfNull(Matrix, $"{nameof(Matrix)} cannot be null, set the matrix");
// //             ArgumentNullException.ThrowIfNull(RightPart, $"{nameof(RightPart)} cannot be null, set the vector");
// //
// //             double alpha, beta;
// //             double squareNorm;
// //
// //             Solution = new double[RightPart.Length];
// //
// //             double[] gglnew = new double[Matrix.Ggl.Length];
// //             double[] ggunew = new double[Matrix.Ggu.Length];
// //             double[] dinew = new double[Matrix.Di.Length];
// //
// //             Matrix.Ggl.Copy(gglnew);
// //             Matrix.Ggu.Copy(ggunew);
// //             Matrix.Di.Copy(dinew);
// //
// //             Vector r = new(RightPart.Length);
// //             Vector z = new(RightPart.Length);
// //             Vector p = new(RightPart.Length);
// //             Vector tmp = new(RightPart.Length);
// //
// //             LU(gglnew, ggunew, dinew);
// //
// //             SparseMatrix.Dot(Matrix, z, p);
// //             
// //             r = Direct(RightPart - (Matrix * Solution), gglnew, dinew);
// //             z = Reverse(r, ggunew);
// //             p = Direct(Matrix * z, gglnew, dinew);
// //
// //             squareNorm = r * r;
// //
// //             for (int iter = 0; iter < MaxIters && squareNorm > Eps; iter++)
// //             {
// //                 alpha = p * r / (p * p);
// //                 squareNorm = (r * r) - (alpha * alpha * (p * p));
// //                 Solution += alpha * z;
// //                 r -= alpha * p;
// //
// //                 tmp = Direct(Matrix * Reverse(r, ggunew), gglnew, dinew);
// //
// //                 beta = -(p * tmp) / (p * p);
// //                 z = Reverse(r, ggunew) + (beta * z);
// //                 p = tmp + (beta * p);
// //             }
// //         }
// //         catch (Exception ex)
// //         {
// //             Console.WriteLine($"We had problem: {ex.Message}");
// //         }
// //     }
// // }
//
// public class LOSLU(int maxIters, double eps) : IterativeSolver(maxIters, eps)
// {
//     public override void Compute()
//     {
//         try
//         {
//             ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
//             ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");
//
//             _solution = new(_vector.Length);
//
//             double[] gglnew = new double[_matrix.Ggl.Length];
//             double[] ggunew = new double[_matrix.Ggu.Length];
//             double[] dinew = new double[_matrix.Di.Length];
//
//             _matrix.Ggl.Copy(gglnew);
//             _matrix.Ggu.Copy(ggunew);
//             _matrix.Di.Copy(dinew);
//
//             Stopwatch sw = Stopwatch.StartNew();
//
//             LU(gglnew, ggunew, dinew);
//
//             var r = Direct(_vector - (_matrix * _solution), gglnew, dinew);
//             var z = Reverse(r, ggunew);
//             var p = Direct(_matrix * z, gglnew, dinew);
//
//             var squareNorm = r * r;
//
//             for (int iter = 0; iter < MaxIters && squareNorm > Eps; iter++)
//             {
//                 var alpha = p * r / (p * p);
//                 squareNorm = (r * r) - (alpha * alpha * (p * p));
//                 _solution += alpha * z;
//                 r -= alpha * p;
//
//                 var tmp = Direct(_matrix * Reverse(r, ggunew), gglnew, dinew);
//
//                 var beta = -(p * tmp) / (p * p);
//                 z = Reverse(r, ggunew) + (beta * z);
//                 p = tmp + (beta * p);
//             }
//
//             sw.Stop();
//
//             _runningTime = sw.Elapsed;
//         }
//         catch (ArgumentNullException ex)
//         {
//             Console.WriteLine($"We had problem: {ex.Message}");
//             throw;
//         }
//         catch (Exception ex)
//         {
//             Console.WriteLine($"We had problem: {ex.Message}");
//         }
//     }
// }
//
//
// // public class CGMCholesky(int maxIterations, double eps) 
// //     : IterativeSolver(maxIterations, eps)
// // {
// //     public override void Compute()
// //     {
// //         ArgumentNullException.ThrowIfNull(Matrix, $"{nameof(Matrix)} cannot be null, set the matrix");
// //         ArgumentNullException.ThrowIfNull(RightPart, $"{nameof(Vector)} cannot be null, set the vector");
// //
// //         double rightPartNorm = RightPart.Norm();
// //         Solution = new double[RightPart.Length];
// //         Array.Fill(Solution, 1.0);
// //         var r = new double[RightPart.Length];
// //         var z = new double[RightPart.Length];
// //         var result = new double[RightPart.Length];
// //         var product = new double[RightPart.Length];
// //         double[] diNew = new double[Matrix.Size];
// //         double[] ggNew = new double[Matrix.Ggl.Length];
// //         Array.Copy(Matrix.Di, diNew, Matrix.Size);
// //         Array.Copy(Matrix.Ggl, ggNew, Matrix.Ggl.Length);
// //         
// //         CholeskyDecomposition(diNew, ggNew);
// //         
// //         SparseMatrix.Dot(Matrix, Solution, product);
// //
// //         for (int i = 0; i < r.Length; i++)
// //         {
// //             r[i] = RightPart[i] - product[i];
// //         }
// //         
// //         MoveForCholesky(diNew, ggNew, r, z);
// //
// //         for (IterationsCount = 0; IterationsCount < MaxIterations && r.Norm() / rightPartNorm >= Eps; IterationsCount++)
// //         {
// //             MoveForCholesky(diNew, ggNew, r, result);
// //             
// //             SparseMatrix.Dot(Matrix, z, product);
// //             double mrDotR = Vector.Dot(result, r);
// //             double alpha = mrDotR / Vector.Dot(product, z);
// //
// //             for (int i = 0; i < Solution.Length; i++)
// //             {
// //                 Solution[i] += alpha * z[i];
// //                 r[i] -= alpha * product[i];
// //             }
// //             
// //             MoveForCholesky(diNew, ggNew, r, result);
// //
// //             double beta = Vector.Dot(result, r) / mrDotR;
// //
// //             for (int i = 0; i < z.Length; i++)
// //             {
// //                 z[i] = result[i] + beta * z[i];
// //             }
// //         }
// //     }
// // }

using System.Collections.Immutable;
using System.Diagnostics;

namespace Boiling.MathHelper;

public static class IterativeSolverExtensions
{
    public static IterativeSolver SetMatrixEx(this IterativeSolver solver, SparseMatrix matrix)
    {
        solver.SetMatrix(matrix);
        return solver;
    }

    public static IterativeSolver SetVectorEx(this IterativeSolver solver, double[] vector)
    {
        solver.SetVector(vector);
        return solver;
    }
}

public abstract class IterativeSolver(int maxIterations, double eps)
{
    protected TimeSpan? _runningTime;
    protected SparseMatrix _matrix = default!;
    protected double[] _vector = default!;
    protected double[]? _solution;

    public int MaxIterations { get; } = maxIterations;
    public double Eps { get; } = eps;
    public TimeSpan? RunningTime => _runningTime;
    public ImmutableArray<double>? Solution => _solution?.ToImmutableArray();

    public void SetMatrix(SparseMatrix matrix) => _matrix = matrix;

    public void SetVector(double[] vector) => _vector = vector;

    public abstract void Compute(double[] initial);

    protected double[] Direct(double[] vector, double[] gglnew, double[] dinew)
    {
        var y = new double[vector.Length];
        Array.Copy(vector, y, vector.Length);
        
        double sum = 0.0;

        for (int i = 0; i < _matrix.Size; i++)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
                sum += gglnew[k] * y[_matrix.Jg[k]];

            y[i] = (y[i] - sum) / dinew[i];
            sum = 0.0;
        }

        return y;
    }

    protected double[] Reverse(double[] vector, double[] ggunew)
    {
        var result = new double[vector.Length];
        Array.Copy(vector, result, vector.Length);

        for (int i = _matrix.Size - 1; i >= 0; i--)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
                result[_matrix.Jg[k]] -= ggunew[k] * result[i];
        }

        return result;
    }

    protected void LU(double[] gglnew, double[] ggunew, double[] dinew)
    {
        double suml = 0.0;
        double sumu = 0.0;
        double sumdi = 0.0;

        for (int i = 0; i < _matrix.Size; i++)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
            {
                int j = _matrix.Jg[k];
                int j0 = _matrix.Ig[j];
                int j1 = _matrix.Ig[j + 1];
                int ik = i0;
                int kj = j0;

                while (ik < k && kj < j1)
                {
                    if (_matrix.Jg[ik] == _matrix.Jg[kj])
                    {
                        suml += gglnew[ik] * ggunew[kj];
                        sumu += ggunew[ik] * gglnew[kj];
                        ik++;
                        kj++;
                    }
                    else if (_matrix.Jg[ik] > _matrix.Jg[kj])
                    {
                        kj++;
                    }
                    else
                    {
                        ik++;
                    }
                }

                gglnew[k] -= suml;
                ggunew[k] = (ggunew[k] - sumu) / dinew[j];
                sumdi += gglnew[k] * ggunew[k];
                suml = 0.0;
                sumu = 0.0;
            }

            dinew[i] -= sumdi;
            sumdi = 0.0;
        }
    }
}

public class LOSLU(int maxIterations, double eps) : IterativeSolver(maxIterations, eps)
{
    public override void Compute(double[] initial)
    {
        try
        {
            ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
            ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");

            _solution = initial;

            double[] gglnew = new double[_matrix.Ggl.Length];
            double[] ggunew = new double[_matrix.Ggu.Length];
            double[] dinew = new double[_matrix.Di.Length];

            Array.Copy(_matrix.Ggl, gglnew, _matrix.Ggl.Length);
            Array.Copy(_matrix.Ggu, ggunew, _matrix.Ggu.Length);
            Array.Copy(_matrix.Di, dinew, _matrix.Di.Length);

            LU(gglnew, ggunew, dinew);

            var tmp = new double[_solution.Length];
            var dotRes = new double[_solution.Length];
            SparseMatrix.Dot(_matrix, _solution, tmp);

            for (int i = 0; i < tmp.Length; i++)
                tmp[i] = _vector[i] - tmp[i];
            
            var r = Direct(tmp, gglnew, dinew);
            var z = Reverse(r, ggunew);
            
            SparseMatrix.Dot(_matrix, z, tmp);
            
            var p = Direct(tmp, gglnew, dinew);

            var squareNorm = r.Norm();
            squareNorm *= squareNorm;

            for (int iter = 0; iter < MaxIterations && squareNorm > Eps; iter++)
            {
                var alpha = p.Dot(r) / p.Dot(p);
                squareNorm = r.Dot(r) - alpha * alpha * p.Dot(p);

                for (int i = 0; i < _solution.Length; i++)
                    _solution[i] += alpha * z[i];

                for (int i = 0; i < r.Length; i++)
                    r[i] -= alpha * p[i];

                tmp = Reverse(r, ggunew);
                SparseMatrix.Dot(_matrix, tmp, dotRes);
                tmp = Direct(dotRes, gglnew, dinew);

                var beta = -p.Dot(tmp) / p.Dot(p);

                for (int i = 0; i < z.Length; i++)
                    dotRes[i] = beta * z[i];
                
                z = Reverse(r, ggunew);
                for (int i = 0; i < dotRes.Length; i++)
                    z[i] += dotRes[i];

                for (int i = 0; i < p.Length; i++)
                    p[i] = tmp[i] + beta * p[i];
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