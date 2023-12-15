using Boiling.Algorithms;
using Boiling.FemContext.BasisInfo;
using Boiling.FemContext.BasisInfo.Interfaces;
using Boiling.MathHelper;
using Boiling.Meshing;

namespace Boiling.FemContext.SLAEAssembler;

public abstract class BaseAssembler
{
    protected readonly SparseMatrix Matrix;
    protected readonly double[] Vector;
    protected readonly Mesh Mesh;
    protected readonly IBasis Basis;
    protected readonly BasisInfoCollection BasisInfo;
    protected readonly TimeMesh TimeMesh;

    public double[] PrevSolution { get; set; } = null!;

    protected BaseAssembler(Mesh mesh, IBasis basis, BasisInfoCollection basisInfo, TimeMesh timeMesh)
    {
        Mesh = mesh;
        Basis = basis;
        BasisInfo = basisInfo;
        TimeMesh = timeMesh;

        var portrait = PortraitBuilder.GeneratePortrait(basisInfo);
        
        Matrix = new SparseMatrix(portrait.Ig, portrait.Jg);
        Vector = new double[portrait.Ig.Length - 1];
    }

    protected abstract void AssembleLocalSlae(int ielem);
    public abstract SparseMatrix GetMatrix(int timeMoment);
    public abstract double[] GetVector(int timeMoment);
}