using System.Diagnostics.CodeAnalysis;
using Boiling.FemContext.BasisInfo;
using Boiling.MathHelper;
using Boiling.Meshing;

namespace Boiling.FemContext.BoundariesHandler;

public class BoundaryHandler(Mesh mesh, BasisInfoCollection basisInfo) 
    : BaseBoundaryHandler(mesh, basisInfo)
{
    private readonly Matrix _localMassMatrix = new(2, 2)
    {
        [0, 0] = 2.0 / 6.0,
        [0, 1] = 1.0 / 6.0,
        [1, 0] = 1.0 / 6.0,
        [1, 1] = 2.0 / 6.0,
    };
    private readonly double[] _thetaVector = new double[2];
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
            
            _thetaVector[0] = n.Theta(p1.X, p1.Y);
            _thetaVector[1] = n.Theta(p2.X, p2.Y);

            double len = Math.Sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y));
            
            Matrix.Dot(_localMassMatrix, _thetaVector, _localVector);

            vector[node1] += len * _localVector[0];
            vector[node2] += len * _localVector[1];
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
            var ubeta1 = n.Ubeta(p1.X, p1.Y);
            var ubeta2 = n.Ubeta(p2.X, p2.Y);
            double len = Math.Sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y));
            
            matrix.Add(node1, node1, beta * len * _localMassMatrix[0, 0]);
            matrix.Add(node1, node2, beta * len * _localMassMatrix[0, 1]);
            matrix.Add(node2, node1, beta * len * _localMassMatrix[1, 0]);
            matrix.Add(node2, node2, beta * len * _localMassMatrix[1, 1]);

            vector[node1] += beta * len * (ubeta1 * _localMassMatrix[0, 0] + ubeta2 * _localMassMatrix[0, 1]);
            vector[node2] += beta * len * (ubeta1 * _localMassMatrix[1, 0] + ubeta2 * _localMassMatrix[1, 1]);
        }
    }
}