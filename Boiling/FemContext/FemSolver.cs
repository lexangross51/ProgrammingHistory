using Boiling.Algorithms;
using Boiling.FemContext.BasisInfo.Interfaces;
using Boiling.FemContext.BoundariesHandler;
using Boiling.FemContext.SLAEAssembler;
using Boiling.MathHelper;
using Boiling.Meshing;

namespace Boiling.FemContext;

public class FemSolver
{
    private readonly BaseAssembler _slaeAssembler;
    private readonly BaseBoundaryHandler _boundaryHandler;
    private readonly DirectSolver _solver;
    private readonly Mesh _mesh;
    private readonly TimeMesh _timeMesh;

    public FemSolver(Mesh mesh, IBasis basis, TimeMesh timeMesh)
    {
        _timeMesh = timeMesh;
        _mesh = mesh;
        
        var basisInfo = Numerator.NumerateBasisFunctions(mesh, basis);

        _slaeAssembler = new Assembler(mesh, basis, basisInfo, timeMesh);
        _boundaryHandler = new BoundaryHandler(mesh, basisInfo);
        _solver = new LUSolver();
    }
    
    private void SetInitialCondition(double temperature)
    {
        _slaeAssembler.PrevSolution = new double[_mesh.Points.Length].Select(_ => temperature).ToArray();
    }
    
    public void Solve()
    {
        SetInitialCondition(25.0);

        var matrix = _slaeAssembler.GetMatrix(1);
        _solver.SetMatrix(matrix);

        for (int i = 1; i < _timeMesh.TimesCount; i++)
        {
            var vector = _slaeAssembler.GetVector(i);

            ApplyBoundaries(matrix, vector, i == 1);
            
            _solver.SetVector(vector);
            _solver.Compute();
            
            Array.Copy(_solver.Solution, _slaeAssembler.PrevSolution, _solver.Solution.Length);
            Console.WriteLine($"Time moment #{i}: {_timeMesh[i]} seconds");
            
            if (i % 100 == 0 || i == _timeMesh.TimesCount - 1)
                Utilities.PrintAtTime(_mesh, _solver.Solution, i, @"C:\Users\lexan\source\repos\Python\meshes");
        }
    }

    private void ApplyBoundaries(SparseMatrix matrix, double[] vector, bool matrixAndVector)
    {
        if (_mesh.Neumann is not null) 
            _boundaryHandler.ApplyNeumann(_mesh.Neumann, vector);
        
        if (_mesh.Newton is not null)
            _boundaryHandler.ApplyNewton(_mesh.Newton, matrix, vector, matrixAndVector);
        
        if (_mesh.Dirichlet.Length != 0)
            _boundaryHandler.ApplyDirichlet(_mesh.Dirichlet, matrix, vector);
    }
}