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
    private readonly IterativeSolver _solver;
    private readonly Mesh _mesh;
    private readonly IBasis _basis;
    // private List<Dirichlet>? _dirichlet;
    // private readonly Newton _newton;

    public IEnumerable<double>? Solution => _solver.Solution;
    
    public FemSolver(Mesh mesh, IBasis basis)
    {
        _basis = basis;
        var basisInfo = Numerator.NumerateBasisFunctions(mesh, basis);

        // Utilities.SaveBasisInfo(mesh, _basisInfo, @"C:\Users\lexan\source\repos\Python");
        _mesh = mesh;
        _slaeAssembler = new Assembler(mesh, basis, basisInfo);
        _boundaryHandler = new BoundaryHandler(mesh, basisInfo);
        _solver = new LOS(10_000, 1E-20);
        
        // RenumerateDirichletNodes();
        // RenumerateNeumannNodes();

        // _newton = new Newton(_mesh, _basis, _basisInfo);

        // Utilities.SaveBiQuadDirichlet(mesh, _basisInfo, _dirichlet!, @"C:\Users\lexan\source\repos\Python");
        // if (_mesh.Neumann != null)
        //     Utilities.SaveBiQuadNeumann(mesh, _basisInfo, _mesh.Neumann, @"C:\Users\lexan\source\repos\Python");
    }

    // private void RenumerateDirichletNodes()
    // {
    //     _dirichlet = new List<Dirichlet>();
    //     
    //     // Edges on the outer border (only dirichlet edges)
    //     var dEdges = _mesh.Elements
    //         .SelectMany(elem => elem.Edges)
    //         .GroupBy(e => e)
    //         .Where(g => g.Count() == 1)
    //         .SelectMany(g => g)
    //         .Where(e => _mesh.Dirichlet.Any(d => d.Node == e.Node1)
    //                     && _mesh.Dirichlet.Any(d => d.Node == e.Node2));
    //
    //     var edgesPerElems = new Dictionary<int, List<Edge>>();
    //     
    //     // Find elements for edges
    //     foreach (var edge in dEdges)
    //     {
    //         int ielem = Array.FindIndex(_mesh.Elements, el => el.Edges.Contains(edge));
    //         
    //         if (!edgesPerElems.ContainsKey(ielem))
    //             edgesPerElems.Add(ielem, new List<Edge>());
    //         
    //         edgesPerElems[ielem].Add(edge);
    //     }
    //     
    //     foreach (var pair in edgesPerElems)
    //     {
    //         int ielem = pair.Key;
    //         var element = _mesh.Elements[ielem];
    //         var edges = pair.Value;
    //
    //         foreach (var edge in edges)
    //         {
    //             int iedge = element.Edges.IndexOf(edge);
    //
    //             var bf = _basisInfo.GetFunctionAtNode(ielem, edge.Node1);
    //             var index = Array.FindIndex(_mesh.Dirichlet, d => d.Node == edge.Node1);
    //             _dirichlet.Add(new Dirichlet(bf.FunctionNumber, _mesh.Dirichlet[index].Value));
    //
    //             bf = _basisInfo.GetFunctionAtNode(ielem, edge.Node2);
    //             index = Array.FindIndex(_mesh.Dirichlet, d => d.Node == edge.Node2);
    //             _dirichlet.Add(new Dirichlet(bf.FunctionNumber, _mesh.Dirichlet[index].Value));
    //             
    //             bf = _basisInfo.GetFunctionAtEdge(ielem, iedge);
    //             _dirichlet.Add(new Dirichlet(bf.FunctionNumber, _mesh.Dirichlet[index].Value));
    //         }
    //     }
    //
    //     _dirichlet = _dirichlet.DistinctBy(d => d.Node).ToList();
    // }
    //
    // private void RenumerateNeumannNodes()
    // {
    //     if (_mesh.Neumann is null) return;
    //     
    //     var edgesPerElems = new Dictionary<int, List<Edge>>();
    //     
    //     // Find elements for edges
    //     foreach (var n in _mesh.Neumann)
    //     {
    //         var edge = n.Border;
    //         int ielem = Array.FindIndex(_mesh.Elements, el => el.Edges.Contains(edge));
    //         
    //         if (!edgesPerElems.ContainsKey(ielem))
    //             edgesPerElems.Add(ielem, new List<Edge>());
    //         
    //         edgesPerElems[ielem].Add(edge);
    //     }
    //     
    //     foreach (var (ielem, edges) in edgesPerElems)
    //     {
    //         foreach (var edge in edges)
    //         {
    //             var bf1 = _basisInfo.GetFunctionAtNode(ielem, edge.Node1);
    //             var bf2 = _basisInfo.GetFunctionAtNode(ielem, edge.Node2);
    //             var index = Array.FindIndex(_mesh.Neumann, d => d.Border.Equal(edge));
    //             var border = _mesh.Neumann[index].Border;
    //             
    //             border.Node1 = bf1.FunctionNumber;
    //             border.Node2 = bf2.FunctionNumber;
    //             
    //             _mesh.Neumann[index].Border = border;
    //         }
    //     }
    // }
    
    public void Solve()
    {
        var slae = _slaeAssembler.GetSlae();

        ApplyBoundaries(slae.Matrix, slae.Vector);
        
        _solver.SetSystem(slae.Matrix, slae.Vector);
        _solver.Compute();
    }

    private void ApplyBoundaries(SparseMatrix matrix, double[] vector)
    {
        if (_mesh.Newton is not null)
            _boundaryHandler.ApplyNewton(_mesh.Newton, matrix, vector);
        
        if (_mesh.Neumann is not null) 
            _boundaryHandler.ApplyNeumann(_mesh.Neumann, vector);
        
        if (_mesh.Dirichlet.Length != 0)
            _boundaryHandler.ApplyDirichlet(_mesh.Dirichlet, matrix, vector);
    }

    // Only for analytical functions
    // [SuppressMessage("ReSharper", "PossibleMultipleEnumeration")]
    // public double RootMeanSquare(IEnumerable<Point> points)
    // {
    //     var func = _dirichlet![0].Value;
    //     double dif = 0.0;
    //     
    //     foreach (var p in points)
    //     {
    //         double exact = func(p.X, p.Y);
    //         // double numeric = ValueAtPoint(p.X, p.Y);
    //
    //         dif += (exact - numeric) * (exact - numeric);
    //     }
    //
    //     return Math.Sqrt(dif / points.Count());
    // }
    
    // public double ValueAtPoint(double x, double y)
    // {
    //     var point = new Point(x, y);
    //     int ielem = FindNumberElement(point);
    //
    //     if (ielem == -1) return double.MinValue;
    //     
    //     var result = 0.0;
    //
    //     // _newton.Point = point;
    //     // _newton.NumberElement = ielem;
    //     //
    //     // _newton.Compute();
    //
    //     for (int i = 0; i < _basis.BasisSize; i++)
    //     {
    //         var p = _newton.Result;
    //         var index = _basisInfo[ielem, i].FunctionNumber;
    //         result += _solver.Solution![index] * _basis.Phi(i, p.X, p.Y);
    //     }
    //
    //     return result;
    // }
    
    //  private int FindNumberElement(Point point)
    // {
    //     const double floatEps = 1E-04;
    //     const double eps = 1E-05;
    //
    //     for (int ielem = 0; ielem < _mesh.Elements.Length; ielem++)
    //     {
    //         var element = _mesh.Elements[ielem];
    //
    //         var intersectionCounter = 0;
    //
    //         var edges = new List<(Point, Point)>
    //         {
    //             (_mesh.Points[element.Nodes[0]], _mesh.Points[element.Nodes[1]]),
    //             (_mesh.Points[element.Nodes[1]], _mesh.Points[element.Nodes[3]]),
    //             (_mesh.Points[element.Nodes[3]], _mesh.Points[element.Nodes[2]]),
    //             (_mesh.Points[element.Nodes[2]], _mesh.Points[element.Nodes[0]])
    //         };
    //
    //         foreach (var edge in edges)
    //         {
    //             if (OnEdge(edge)) return ielem;
    //
    //             var pt = point;
    //             (Point A, Point B) sortedEdge = edge.Item1.Y > edge.Item2.Y ? (edge.Item2, edge.Item1) : edge;
    //
    //             if (Math.Abs(pt.Y - sortedEdge.A.Y) < floatEps || Math.Abs(pt.Y - sortedEdge.B.Y) < floatEps)
    //             {
    //                 pt = pt with { Y = pt.Y + eps };
    //             }
    //
    //             var maxX = edge.Item1.X > edge.Item2.X ? edge.Item1.X : edge.Item2.X;
    //             var minX = edge.Item1.X < edge.Item2.X ? edge.Item1.X : edge.Item2.X;
    //
    //             if (pt.Y > sortedEdge.B.Y || pt.Y < sortedEdge.A.Y || pt.X > maxX) continue;
    //
    //             if (pt.X < minX)
    //             {
    //                 intersectionCounter++;
    //             }
    //             else
    //             {
    //                 var mRed = Math.Abs(sortedEdge.A.X - sortedEdge.B.X) > uint.MinValue
    //                     ? (sortedEdge.B.Y - sortedEdge.A.Y) / (sortedEdge.B.X - sortedEdge.A.X)
    //                     : double.PositiveInfinity;
    //                 var mBlue = Math.Abs(sortedEdge.A.X - pt.X) > uint.MinValue
    //                     ? (pt.Y - sortedEdge.A.Y) / (pt.X - sortedEdge.A.X)
    //                     : double.PositiveInfinity;
    //
    //                 if (mBlue >= mRed) intersectionCounter++;
    //             }
    //         }
    //
    //         if (intersectionCounter % 2 == 1) return ielem;
    //
    //         bool OnEdge((Point, Point) edge)
    //         {
    //             var distance1 = Math.Sqrt((point.X - edge.Item1.X) * (point.X - edge.Item1.X) +
    //                                       (point.Y - edge.Item1.Y) * (point.Y - edge.Item1.Y));
    //             var distance2 = Math.Sqrt((point.X - edge.Item2.X) * (point.X - edge.Item2.X) +
    //                                       (point.Y - edge.Item2.Y) * (point.Y - edge.Item2.Y));
    //             var edgeLength = Math.Sqrt((edge.Item2.X - edge.Item1.X) * (edge.Item2.X - edge.Item1.X) +
    //                                        (edge.Item2.Y - edge.Item1.Y) * (edge.Item2.Y - edge.Item1.Y));
    //
    //             return Math.Abs(distance1 + distance2 - edgeLength) < floatEps;
    //         }
    //     }
    //
    //     return -1;
    // }
}