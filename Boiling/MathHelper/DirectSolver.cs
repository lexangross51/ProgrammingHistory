namespace Boiling.MathHelper;

public abstract class DirectSolver
{
    protected ProfileMatrix Matrix = null!;
    protected double[] Vector = null!;
    public double[] Solution { get; private set; } = null!;
    protected abstract void Decompose();
    
    public void SetMatrix(SparseMatrix matrix)
    {
        Matrix = matrix.ToProfile();
        Solution = new double[Matrix.Size];
        Decompose();
    }

    public void SetVector(double[] vector)
    {
        Vector = vector;
    }

    public abstract void Compute();
}

public class LUSolver : DirectSolver
{
    public override void Compute()
    {
        Array.Clear(Solution);
        
        for (int i = 0; i < Matrix.Size; i++)
        {
            int i0 = Matrix.Ig[i];
            int i1 = Matrix.Ig[i + 1];
            int j = i - (i1 - i0);
            double sum = 0.0;

            for (int k = i0; k < i1; k++)
            {
                sum += Matrix.Ggl[k] * Solution[j++];
            }
        
            Solution[i] = (Vector[i] - sum) / Matrix.Di[i];
        }

        for (int i = Matrix.Size - 1; i >= 0; i--)
        {
            int i0 = Matrix.Ig[i];
            int i1 = Matrix.Ig[i + 1];
            int j = i - (i1 - i0);

            for (int k = i0; k < i1; k++)
            {
                Solution[j++] -= Matrix.Ggu[k] * Solution[i];
            }
        }
    }

    protected override void Decompose()
    {
        for (int i = 1; i < Matrix.Size; i++)
        {
            int j0 = i - (Matrix.Ig[i + 1] - Matrix.Ig[i]);
            
            for (int ii = Matrix.Ig[i]; ii < Matrix.Ig[i + 1]; ii++)
            {
                int j = ii - Matrix.Ig[i] + j0;
                double sumL = 0, sumU = 0;
                
                if (Matrix.Ig[j] < Matrix.Ig[j + 1])
                {
                    int j0J = j - (Matrix.Ig[j + 1] - Matrix.Ig[j]);
                    int jjBeg = j0 < j0J ? j0J : j0;
                    int jjEnd = j < i - 1 ? j : i - 1;
                    
                    for (int k = 0; k < jjEnd - jjBeg; k++)
                    {
                        int indPrev = Matrix.Ig[j] + jjBeg - j0J + k;
                        int indNow = Matrix.Ig[i] + jjBeg - j0 + k;
                        sumL += Matrix.Ggu[indPrev] * Matrix.Ggl[indNow];
                        sumU += Matrix.Ggu[indNow] * Matrix.Ggl[indPrev];
                    }
                }
                Matrix.Ggl[ii] -= sumL;
                Matrix.Ggu[ii] -= sumU;
                Matrix.Ggu[ii] /= Matrix.Di[j];
                Matrix.Di[i] -= Matrix.Ggl[ii] * Matrix.Ggu[ii];
            }
        }
    }
}