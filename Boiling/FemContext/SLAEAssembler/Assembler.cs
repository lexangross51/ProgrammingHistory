using Boiling.FemContext.BasisInfo;
using Boiling.FemContext.BasisInfo.Interfaces;
using Boiling.MathHelper;
using Boiling.Meshing;
using Boiling.Meshing.Geometry;

namespace Boiling.FemContext.SLAEAssembler;

public class Assembler(Mesh mesh, IBasis basis, BasisInfoCollection basisInfo, TimeMesh timeMesh) 
    : BaseAssembler(mesh, basis, basisInfo, timeMesh)
{
    private const double Velocity = 0.001;
    private readonly Matrix _stiffnessMatrix = new(basis.BasisSize, basis.BasisSize);
    private readonly Matrix _massMatrix = new(basis.BasisSize, basis.BasisSize);
    private readonly Matrix _navierStokesMatrix = new(basis.BasisSize, basis.BasisSize);
    private readonly double[] _localVector = new double[basis.BasisSize];
    private readonly double[] _localRightPart = new double[basis.BasisSize];
    // private readonly Integration _integrator = new(Quadratures.GaussOrder3());
    // private readonly Vector<double> _gradPhiI = new(2);
    // private readonly Vector<double> _gradPhiJ = new(2);
    // private readonly double[] _matrixGradI = new double[2];
    // private readonly double[] _matrixGradJ = new double[2];
    // private readonly Rectangle _masterElement = new(new Point(), new Point(1, 1));
    private readonly double _radius = mesh.Points.Max(p => p.X);
    private readonly double _height = mesh.Points.Max(p => p.Y);

    private Point CalculateVelocity(int ielem)
    {
        //    2
        // 3     1
        //    4
        var element = Mesh.Elements[ielem]; 
        var nodes = element.Nodes;
        double x = (Mesh.Points[nodes[0]].X + Mesh.Points[nodes[1]].X) * 0.5; 
        double y = (Mesh.Points[nodes[0]].Y + Mesh.Points[nodes[2]].Y) * 0.5; 

        return element.VelocityArea switch
        {
            1 => x > 0.75 * _radius
                ? new Point(0.0, Velocity * (_radius - x) / (_radius * 0.25))
                : new Point(0.0, Velocity * (x - _radius * 0.5) / (_radius * 0.25)),
            2 => y > 0.75 * _height
                ? new Point(-Velocity * (_height - y) / (_height * 0.25), 0.0)
                : new Point(-Velocity * (y - _height * 0.5) / (_height * 0.25), 0.0),
            3 => x < _radius * 0.25
                ? new Point(0.0, -Velocity * x / (_radius * 0.25))
                : new Point(0.0, -Velocity * (_radius * 0.5 - x) / (_radius * 0.25)),
            4 => y < _height * 0.25
                ? new Point(Velocity * y / (_height * 0.25), 0.0)
                : new Point(Velocity * (_height * 0.5 - y) / (_height * 0.25), 0.0),
            _ => new Point()
        };
    }
    
    protected override void AssembleLocalSlae(int ielem)
    {
        var element = Mesh.Elements[ielem];
        var nodes = element.Nodes;
        var r = Mesh.Points[nodes[0]].X;
        var r1 = Mesh.Points[nodes[^1]].X;
        var z = Mesh.Points[nodes[0]].Y;
        var z1 = Mesh.Points[nodes[^1]].Y;
        var hr = r1 - r;
        var hz = z1 - z;
        var v = CalculateVelocity(ielem);
        
        // Assembly local stiffness matrix
        #region Numerical integration

        // for (int i = 0; i < Basis.BasisSize; i++)
        // {
        //     int i1 = i;
        //     
        //     for (int j = 0; j <= i; j++)
        //     {
        //         int j1 = j;
        //
        //         double ScalarFunc(double ksi, double eta)
        //         {
        //             var jacobiMatrix = FemHelper.CalculateJacobiMatrix(Mesh, ielem, Basis, BasisInfo, ksi, eta);
        //             double jacobian = FemHelper.Jacobian(jacobiMatrix);
        //             FemHelper.InvertJacobiMatrix(jacobiMatrix);
        //             
        //             _gradPhiI[0] = Basis.DPhi(i1, 0, ksi, eta);
        //             _gradPhiI[1] = Basis.DPhi(i1, 1, ksi, eta);
        //             _gradPhiJ[0] = Basis.DPhi(j1, 0, ksi, eta);
        //             _gradPhiJ[1] = Basis.DPhi(j1, 1, ksi, eta);
        //
        //             MathHelper.Matrix.Dot(jacobiMatrix, _gradPhiI, _matrixGradI);
        //             MathHelper.Matrix.Dot(jacobiMatrix, _gradPhiJ, _matrixGradJ);
        //             
        //             return (_matrixGradI[0] * _matrixGradJ[0] + _matrixGradI[1] * _matrixGradJ[1]) *
        //                    Math.Abs(jacobian) * (r + hr * ksi);
        //         }
        //         
        //         _stiffnessMatrix[i, j] = _stiffnessMatrix[j, i] = _integrator.Integrate2D(ScalarFunc, _masterElement);
        //     }
        // }

        #endregion

        #region Analytical

        double g1 = hz / (12.0 * hr) * (r + r1);
        double g2 = hr / (12.0 * hz) * (3.0 * r + r1);
        double g3 = hr / (12.0 * hz) * (r + 3.0 * r1);
        double g4 = hr / (12.0 * hz) * (r + r1);
        
        _stiffnessMatrix[0, 0] = 2.0 * g1 + g2;
        _stiffnessMatrix[0, 1] = _stiffnessMatrix[1, 0] = -2.0 * g1 + g4;
        _stiffnessMatrix[0, 2] = _stiffnessMatrix[2, 0] = g1 - g2;
        _stiffnessMatrix[0, 3] = _stiffnessMatrix[3, 0] = -g1 - g4;
        _stiffnessMatrix[1, 1] = 2.0 * g1 + g3;
        _stiffnessMatrix[1, 2] = _stiffnessMatrix[2, 1] = -g1 - g4;
        _stiffnessMatrix[1, 3] = _stiffnessMatrix[3, 1] = g1 - g3;
        _stiffnessMatrix[2, 2] = 2.0 * g1 + g2;
        _stiffnessMatrix[2, 3] = _stiffnessMatrix[3, 2] = -2.0 * g1 + g4;
        _stiffnessMatrix[3, 3] = 2.0 * g1 + g3;
        
        #endregion
        
        // Assembly local mass matrix
        #region Numerical integration

        // for (int i = 0; i < Basis.BasisSize; i++)
        // {
        //     int i1 = i;
        //
        //     for (int j = 0; j <= i; j++)
        //     {
        //         int j1 = j;
        //
        //         double ScalarFunc(double ksi, double eta)
        //         {
        //             var jacobiMatrix = FemHelper.CalculateJacobiMatrix(Mesh, ielem, Basis, BasisInfo, ksi, eta);
        //             double jacobian = FemHelper.Jacobian(jacobiMatrix);
        //
        //             return Basis.Phi(i1, ksi, eta) * Basis.Phi(j1, ksi, eta) * Math.Abs(jacobian) * (r + hr * ksi);
        //         }
        //
        //         _massMatrix[i, j] = _massMatrix[j, i] = _integrator.Integrate2D(ScalarFunc, _masterElement);
        //     }
        // }

        #endregion

        #region Analytical

        double m1 = hr * hz * (3.0 * r + r1) / 72.0;
        double m2 = hr * hz * (r + r1) / 72.0;
        double m3 = hr * hz * (r + 3.0 * r1) / 72.0;
        
        _massMatrix[0, 0] = 2.0 * m1;
        _massMatrix[0, 1] = _massMatrix[1, 0] = 2.0 * m2;
        _massMatrix[0, 2] = _massMatrix[2, 0] = m1;
        _massMatrix[0, 3] = _massMatrix[3, 0] = m2;
        _massMatrix[1, 1] = 2.0 * m3;
        _massMatrix[1, 2] = _massMatrix[2, 1] = m2;
        _massMatrix[1, 3] = _massMatrix[3, 1] = m3;
        _massMatrix[2, 2] = 2.0 * m1;
        _massMatrix[2, 3] = _massMatrix[3, 2] = 2.0 * m2;
        _massMatrix[3, 3] = 2.0 * m3;

        #endregion
        
        // Assembly local right part vector
        var source = Mesh.Materials[Mesh.Elements[ielem].AreaNumber].Source;
        for (int i = 0; i < Basis.BasisSize; i++)
        {
            var bf = BasisInfo[ielem, i];
            var x = Mesh.Points[bf.Index].X;
            var y = Mesh.Points[bf.Index].Y;
        
            _localRightPart[i] = source(x, y);
        }
        
        for (int i = 0; i < Basis.BasisSize; i++)
        {
            _localVector[i] = 0.0;
            
            for (int j = 0; j < Basis.BasisSize; j++)
            {
                _localVector[i] += _massMatrix[i, j] * _localRightPart[j];
            }
        }

        // Assembly local navier stokes matrix
        #region Numerical integration

        // for (int i = 0; i < Basis.BasisSize; i++)
        // {
        //     var i1 = i;
        //     
        //     for (int j = 0; j < Basis.BasisSize; j++)
        //     {
        //         var j1 = j;
        //
        //         double ScalarFunc(double ksi, double eta)
        //         {
        //             var jacobiMatrix = FemHelper.CalculateJacobiMatrix(Mesh, ielem, Basis, BasisInfo, ksi, eta);
        //             var jacobian = FemHelper.Jacobian(jacobiMatrix);
        //             FemHelper.InvertJacobiMatrix(jacobiMatrix);
        //             
        //             _gradPhiI[0] = Basis.DPhi(j1, 0, ksi, eta);
        //             _gradPhiI[1] = Basis.DPhi(j1, 1, ksi, eta);
        //             
        //             MathHelper.Matrix.Dot(jacobiMatrix, _gradPhiI, _matrixGradI);
        //
        //             return (_matrixGradI[0] * v.X + _matrixGradI[1] * v.Y) * Basis.Phi(i1, ksi, eta) *
        //                    Math.Abs(jacobian) * (r + hr * ksi);
        //         }
        //
        //         _navierStokesMatrix[i, j] = _integrator.Integrate2D(ScalarFunc, _masterElement);
        //     }
        // }

        #endregion

        #region Analytical

        _navierStokesMatrix[0, 0] = r * r * v.Y * 0.125 + (-6.0 * r1 * v.Y - 8.0 * v.X * hz) * r / 72.0 - (r1 * v.Y + 4.0 / 3.0 * v.X * hz) * r1 / 24.0;
        _navierStokesMatrix[0, 1] = hz * (r + 0.5 * r1) * v.X / 9.0 + r * r * v.Y / 24.0 - r1 * r1 * v.Y / 24.0;
        _navierStokesMatrix[0, 2] = -r * r * v.Y * 0.125 + (6.0 * r1 * v.Y - 4.0 * v.X * hz) * r / 72.0 + (r1 * v.Y - 2.0 / 3.0 * v.X * hz) * r1 / 24.0;
        _navierStokesMatrix[0, 3] = hz * (r + 0.5 * r1) * v.X / 18.0 - r * r * v.Y / 24.0 + r1 * r1 * v.Y / 24.0;

        _navierStokesMatrix[1, 0] = -hz * (r + 2.0 * r1) * v.X / 18.0 + r * r * v.Y / 24.0 - r1 * r1 * v.Y / 24.0;
        _navierStokesMatrix[1, 1] = r * r * v.Y / 24.0 + (6.0 * r1 * v.Y + 4.0 * v.X * hz) * r / 72.0 - (r1 * v.Y - 8.0 / 9.0 * v.X * hz) * r1 * 0.125;
        _navierStokesMatrix[1, 2] = -hz * (r + 2.0 * r1) * v.X / 36.0 - r * r * v.Y / 24.0 + r1 * r1 * v.Y / 24.0;
        _navierStokesMatrix[1, 3] = -r * r * v.Y / 24.0 + (-6.0 * r1 * v.Y + 2.0 * v.X * hz) * r / 72.0 + (r1 * v.Y + 4.0 / 9.0 * v.X * hz) * r1 * 0.125;

        _navierStokesMatrix[2, 0] = r * r * v.Y * 0.125 + (-6.0 * r1 * v.Y - 4.0 * v.X * hz) * r / 72.0 - (r1 * v.Y + 2.0 / 3.0 * v.X * hz) * r1 / 24.0;
        _navierStokesMatrix[2, 1] = hz * (r + 0.5 * r1) * v.X / 18.0 + r * r * v.Y / 24.0 - r1 * r1 * v.Y / 24.0;
        _navierStokesMatrix[2, 2] = -r * r * v.Y * 0.125 + (6.0 * r1 * v.Y - 8.0 * v.X * hz) * r / 72.0 + r1 * (r1 * v.Y - 4.0 / 3.0 * v.X * hz) / 24.0;
        _navierStokesMatrix[2, 3] = hz * (r + 0.5 * r1) * v.X / 9.0 - r * r * v.Y / 24.0 + r1 * r1 * v.Y / 24.0;

        _navierStokesMatrix[3, 0] = -hz * (r + 2.0 * r1) * v.X / 36.0 + r * r * v.Y / 24.0 - r1 * r1 * v.Y / 24.0;
        _navierStokesMatrix[3, 1] = r * r * v.Y / 24.0 + (6.0 * r1 * v.Y + 2.0 * v.X * hz) * r / 72.0 - (r1 * v.Y - 4.0 / 9.0 * v.X * hz) * r1 * 0.125;
        _navierStokesMatrix[3, 2] = -hz * (r + 2.0 * r1) * v.X / 18.0 - r * r * v.Y / 24.0 + r1 * r1 * v.Y / 24.0;
        _navierStokesMatrix[3, 3] = -r * r * v.Y / 24.0 + (-6.0 * r1 * v.Y + 4.0 * v.X * hz) * r / 72.0 + (r1 * v.Y + 8.0 / 9.0 * v.X * hz) * r1 * 0.125;

        #endregion
    }

    public override SparseMatrix GetMatrix(int timeMoment)
    {
        // Utilities.PrintVelocities(Mesh, CalculateVelocity, @"C:\Users\lexan\source\repos\Python\meshes");
        Matrix.Clear();

        double t0 = TimeMesh[timeMoment - 1];
        double t1 = TimeMesh[timeMoment];
        double dt = t1 - t0;

        for (int ielem = 0; ielem < Mesh.Elements.Length; ielem++)
        {
            double lambda = Mesh.Materials[Mesh.Elements[ielem].AreaNumber].Lambda;
            double rhoCp = Mesh.Materials[Mesh.Elements[ielem].AreaNumber].RhoCp;
            
            AssembleLocalSlae(ielem);

            for (int i = 0; i < Basis.BasisSize; i++)
            {
                int globalI = BasisInfo[ielem, i].FunctionNumber;
                    
                for (int j = 0; j < Basis.BasisSize; j++)
                {
                    int globalJ = BasisInfo[ielem, j].FunctionNumber;

                    Matrix.Add(globalI, globalJ, 
                        lambda * _stiffnessMatrix[i, j] + 
                        rhoCp / dt * _massMatrix[i, j] + 
                        rhoCp * _navierStokesMatrix[i, j]);
                }
            }
        }

        return Matrix;
    }

    public override double[] GetVector(int timeMoment)
    {
        Array.Fill(Vector, 0.0);

        double t0 = TimeMesh[timeMoment - 1];
        double t1 = TimeMesh[timeMoment];
        double dt = t1 - t0;
        
        for (int ielem = 0; ielem < Mesh.Elements.Length; ielem++)
        {
            double rhoCp = Mesh.Materials[Mesh.Elements[ielem].AreaNumber].RhoCp;
            
            AssembleLocalSlae(ielem);

            for (int i = 0; i < Basis.BasisSize; i++)
            {
                int globalI = BasisInfo[ielem, i].FunctionNumber;
                Vector[globalI] += _localVector[i];

                for (int j = 0; j < Basis.BasisSize; j++)
                {
                    int globalJ = BasisInfo[ielem, j].FunctionNumber;
                    Vector[globalI] += rhoCp / dt * _massMatrix[i, j] * PrevSolution[globalJ];
                }
            }
        }
        
        return Vector;
    }
}