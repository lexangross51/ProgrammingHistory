using System.Diagnostics.CodeAnalysis;
using Boiling.FemContext.BasisInfo;
using Boiling.MathHelper;
using Boiling.Meshing;

namespace Boiling.FemContext.BoundariesHandler;

public class BoundaryHandler(Mesh mesh, BasisInfoCollection basisInfo) 
    : BaseBoundaryHandler(mesh, basisInfo)
{
    private readonly Matrix[] _localMassMatrixRz = 
    {
        new(2, 2)
        {
            [0, 0] = 2.0 / 6.0,
            [0, 1] = 1.0 / 6.0,
            [1, 0] = 1.0 / 6.0,
            [1, 1] = 2.0 / 6.0,
        },
        new(2, 2)
        {
            [0, 0] = 1.0 / 12.0,
            [0, 1] = 1.0 / 12.0,
            [1, 0] = 1.0 / 12.0,
            [1, 1] = 3.0 / 12.0,
        }
    };
    private readonly Matrix _totalLocalMassMatrix = new(2, 2);
    private readonly double[] _flowVector = new double[2];
    private readonly double[] _localVector = new double[2];

    [SuppressMessage("ReSharper", "PossibleMultipleEnumeration")]
    public override void ApplyDirichlet(IEnumerable<Dirichlet> dirichlet, SparseMatrix matrix, double[] vector)
    {
        foreach (var d in dirichlet)
        {
            var p = Mesh.Points[d.Node];
            matrix.Di[d.Node] = 1.0;
            vector[d.Node] = d.Value(p.X, p.Y);

            int i0 = matrix.Ig[d.Node];
            int i1 = matrix.Ig[d.Node + 1];

            for (int k = i0; k < i1; k++)
                matrix.Ggl[k] = 0.0;

            for (int k = d.Node + 1; k < matrix.Size; k++)
            {
                i0 = matrix.Ig[k];
                i1 = matrix.Ig[k + 1];

                for (int j = i0; j < i1; j++)
                {
                    if (matrix.Jg[j] == d.Node)
                        matrix.Ggu[j] = 0.0;
                }
            }   
        }
    }

    public override void ApplyNeumann(IEnumerable<Neumann> neumann, double[] vector)
    {
        foreach (var n in neumann)
        {
            int node1 = n.Border.Node1;
            int node2 = n.Border.Node2;

            var p1 = Mesh.Points[node1];
            var p2 = Mesh.Points[node2];
            var h = Math.Sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y));
            
            _flowVector[0] = n.Theta(p1.X, p1.Y);
            _flowVector[1] = n.Theta(p2.X, p2.Y);

            _totalLocalMassMatrix[0, 0] = h * (p1.X * _localMassMatrixRz[0][0, 0] + h * _localMassMatrixRz[1][0, 0]);
            _totalLocalMassMatrix[0, 1] = h * (p1.X * _localMassMatrixRz[0][0, 1] + h * _localMassMatrixRz[1][0, 1]);
            _totalLocalMassMatrix[1, 0] = h * (p1.X * _localMassMatrixRz[0][1, 0] + h * _localMassMatrixRz[1][1, 0]);
            _totalLocalMassMatrix[1, 1] = h * (p1.X * _localMassMatrixRz[0][1, 1] + h * _localMassMatrixRz[1][1, 1]);
            
            Matrix.Dot(_totalLocalMassMatrix, _flowVector, _localVector);

            vector[node1] += _localVector[0];
            vector[node2] += _localVector[1];
        }
    }

    public override void ApplyNewton(IEnumerable<Newton> newton, SparseMatrix matrix, double[] vector, 
        bool matrixAndVector = false)
    {
        foreach (var n in newton)
        {
            var node1 = n.Border.Node1;
            var node2 = n.Border.Node2;
            var p1 = Mesh.Points[node1];
            var p2 = Mesh.Points[node2];
            var beta = n.Beta;
            double h = Math.Sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y));
            
            _flowVector[0] = n.Ubeta(p1.X, p1.Y);
            _flowVector[1] = n.Ubeta(p2.X, p2.Y);
            
            _totalLocalMassMatrix[0, 0] = beta * h * (p1.X * _localMassMatrixRz[0][0, 0] + h * _localMassMatrixRz[1][0, 0]);
            _totalLocalMassMatrix[0, 1] = beta * h * (p1.X * _localMassMatrixRz[0][0, 1] + h * _localMassMatrixRz[1][0, 1]);
            _totalLocalMassMatrix[1, 0] = beta * h * (p1.X * _localMassMatrixRz[0][1, 0] + h * _localMassMatrixRz[1][1, 0]);
            _totalLocalMassMatrix[1, 1] = beta * h * (p1.X * _localMassMatrixRz[0][1, 1] + h * _localMassMatrixRz[1][1, 1]);

            Matrix.Dot(_totalLocalMassMatrix, _flowVector, _localVector);

            if (matrixAndVector)
            {
                matrix.Add(node1, node1, _totalLocalMassMatrix[0, 0]);
                matrix.Add(node1, node2, _totalLocalMassMatrix[0, 1]);
                matrix.Add(node2, node1, _totalLocalMassMatrix[1, 0]);
                matrix.Add(node2, node2, _totalLocalMassMatrix[1, 1]);   
            }

            vector[node1] += _localVector[0];
            vector[node2] += _localVector[1];
        }
    }
}