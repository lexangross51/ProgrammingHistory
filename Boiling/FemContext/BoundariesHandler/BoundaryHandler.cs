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
        var fCount = BasisInfo.FunctionsCount;
        var bc1 = new int[fCount];
        
        Array.Fill(bc1, -1);

        var dir = dirichlet.ToArray();
        
        for (int i = 0; i < dir.Length; i++)
        {
            bc1[dir[i].Node] = i;
        }

        for (int i = 0; i < fCount; i++)
        {
            int k;

            if (bc1[i] != -1)
            {
                var (node, value) = dir[bc1[i]];

                matrix.Di[i] = 1.0;

                var p = Mesh.Points[node];
                vector[i] = value(p.X, p.Y);

                for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
                {
                    k = matrix.Jg[j];
                    if (bc1[k] == -1)
                    {
                        vector[k] -= matrix.Ggl[j] * vector[i];
                    }
                    matrix.Ggl[j] = 0.0;
                    matrix.Ggu[j] = 0.0;
                }
            }
            else
            {
                for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
                {
                    k = matrix.Jg[j];
                    if (bc1[k] != -1)
                    {
                        vector[i] -= matrix.Ggl[j] * vector[k];
                        matrix.Ggl[j] = 0.0;
                        matrix.Ggu[j] = 0.0;
                    }
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

    public override void ApplyNewton(IEnumerable<Newton> newton, SparseMatrix matrix, double[] vector)
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

            matrix.Add(node1, node1, _totalLocalMassMatrix[0, 0]);
            matrix.Add(node1, node2, _totalLocalMassMatrix[0, 1]);
            matrix.Add(node2, node1, _totalLocalMassMatrix[1, 0]);
            matrix.Add(node2, node2, _totalLocalMassMatrix[1, 1]);

            vector[node1] += _localVector[0];
            vector[node2] += _localVector[1];
        }
    }
}