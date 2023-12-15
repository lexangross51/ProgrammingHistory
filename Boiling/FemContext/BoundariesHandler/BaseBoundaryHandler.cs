using Boiling.FemContext.BasisInfo;
using Boiling.MathHelper;
using Boiling.Meshing;

namespace Boiling.FemContext.BoundariesHandler;

public abstract class BaseBoundaryHandler(Mesh mesh, BasisInfoCollection basisInfo)
{
    protected readonly Mesh Mesh = mesh;
    protected readonly BasisInfoCollection BasisInfo = basisInfo;

    public abstract void ApplyDirichlet(IEnumerable<Dirichlet> dirichlet, SparseMatrix matrix, double[] vector);
    public abstract void ApplyNeumann(IEnumerable<Neumann> neumann, double[] vector);
    public abstract void ApplyNewton(IEnumerable<Newton> newton, SparseMatrix matrix, double[] vector,
        bool matrixAndVector = false);
}