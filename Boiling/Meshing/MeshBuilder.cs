using Boiling.Meshing.Geometry;
using Boiling.Meshing.Interfaces;

namespace Boiling.Meshing;

public class MeshBuilder : IMeshBuilder
{
    private readonly MeshParameters _meshParameters;
    private readonly Point[] _points;
    private readonly int[] _ix;
    private readonly int[] _iy;
    private readonly FiniteElement[] _elements;
    private List<int> _fictitiousNodes;
    private List<int> _fictitiousElements;
    
    private HashSet<Dirichlet> _dirichlet;
    private HashSet<Neumann>? _neumann;
    private HashSet<Newton>? _newton;

    public IEnumerable<Point> Points { get; set; } = null!;

    public MeshBuilder(MeshParameters meshParameters)
    {
        _meshParameters = meshParameters;
        
        PrepareRefinement();
        
        var totalNx = meshParameters.AbscissaSplits.Sum();
        var totalNy = meshParameters.OrdinateSplits.Sum();
        
        _points = new Point[(totalNx + 1) * (totalNy + 1)];
        _elements = new FiniteElement[totalNx * totalNy];
        _ix = new int[_meshParameters.AbscissaPointsCount];
        _iy = new int[_meshParameters.OrdinatePointsCount];
        _fictitiousNodes = new List<int>();
        _fictitiousElements = new List<int>();

        _dirichlet = new HashSet<Dirichlet>();
    }

    private void PrepareRefinement()
    {
        if (_meshParameters.Refinement == 0) return;

        for (int i = 0; i < _meshParameters.AbscissaSplits.Length; i++)
        {
            _meshParameters.AbscissaSplits[i] *= (int)Math.Pow(2, _meshParameters.Refinement);
            _meshParameters.AbscissaK[i] = Math.Sign(_meshParameters.AbscissaK[i]) * 
                                           Math.Pow(Math.Abs(_meshParameters.AbscissaK[i]),
                                               1.0 / Math.Pow(2, _meshParameters.Refinement));
        }
        
        for (int i = 0; i < _meshParameters.OrdinateSplits.Length; i++)
        {
            _meshParameters.OrdinateSplits[i] *= (int)Math.Pow(2, _meshParameters.Refinement);
            _meshParameters.OrdinateK[i] = Math.Sign(_meshParameters.OrdinateK[i]) * 
                                           Math.Pow(Math.Abs(_meshParameters.OrdinateK[i]),
                                               1.0 / Math.Pow(2, _meshParameters.Refinement));
        }
    }

    public void CreatePoints()
    {
        int totalNx = _meshParameters.AbscissaSplits.Sum() + 1;  // Nodes count on abscissa with splits
        int totalNy = _meshParameters.OrdinateSplits.Sum() + 1;  // Nodes count on ordinate with splits

        var primaryNx = _meshParameters.AbscissaPointsCount;   // Primary abscissa points count
        var primaryNy = _meshParameters.OrdinatePointsCount;   // Primary ordinate points count
        
        int ordinateSplits = 0;
        int abscissaSplits;
        
        // Forming nodes on main horizontal lines
        for (int i = 0; i < primaryNy; i++)
        {
            abscissaSplits = 0;
            
            for (int j = 0; j < primaryNx - 1; j++)
            {
                var a = _meshParameters.ControlPoints[i * primaryNx + j];       // Start
                var b = _meshParameters.ControlPoints[i * primaryNx + j + 1];   // End
                int splits = _meshParameters.AbscissaSplits[j];
                double k = _meshParameters.AbscissaK[j] < 0
                    ? -1.0 / _meshParameters.AbscissaK[j]
                    : _meshParameters.AbscissaK[j];
                double h = (Math.Abs(k - 1.0) < 1E-14)
                    ? (b.X - a.X) / splits
                    : (b.X - a.X) * (1.0 - k) / (1.0 - Math.Pow(k, splits));

                _points[totalNx * ordinateSplits + abscissaSplits] = a;
                
                for (int l = 1; l < splits + 1; l++)
                {
                    double x = _points[totalNx * ordinateSplits + abscissaSplits + l - 1].X + h * Math.Pow(k, l - 1);
                    double t = (x - a.X) / (b.X - a.X);
                    double y = a.Y + t * (b.Y - a.Y);
                    
                    _points[totalNx * ordinateSplits + abscissaSplits + l] = new Point(x, y);
                }

                abscissaSplits += splits;
            }

            if (i < primaryNy - 1)
                ordinateSplits += _meshParameters.OrdinateSplits[i];
        }
        
        // Forming nodes on main vertical lines
        abscissaSplits = 0;
        
        for (int i = 0; i < primaryNx; i++)
        {
            ordinateSplits = 0;
            
            for (int j = 0; j < primaryNy - 1; j++)
            {
                var a = _meshParameters.ControlPoints[j * primaryNx + i];
                var b = _meshParameters.ControlPoints[(j + 1) * primaryNx + i];
                int splits = _meshParameters.OrdinateSplits[j];
                double k = _meshParameters.OrdinateK[j] < 0
                    ? -1.0 / _meshParameters.OrdinateK[j]
                    : _meshParameters.OrdinateK[j];
                double h = Math.Abs(k - 1.0) < 1E-14
                    ? (b.Y - a.Y) / splits
                    : (b.Y - a.Y) * (1.0 - k) / (1.0 - Math.Pow(k, splits));

                _points[totalNx * ordinateSplits + abscissaSplits] = a;
                
                for (int l = 1; l < splits + 1; l++)
                {
                    double y = _points[totalNx * ordinateSplits + abscissaSplits + (l - 1) * totalNx].Y + h * Math.Pow(k, l - 1);
                    double t = (y - a.Y) / (b.Y - a.Y);
                    double x = a.X + t * (b.X - a.X);
                    
                    _points[totalNx * ordinateSplits + abscissaSplits + l * totalNx] = new Point(x, y);
                }

                ordinateSplits += splits;
            }
            
            if (i < primaryNx - 1)
                abscissaSplits += _meshParameters.AbscissaSplits[i];
        }
        
        // Form inner nodes
        for (int i = 1; i < totalNy; i++)
        {
            abscissaSplits = 0;
            
            for (int j = 0; j < primaryNx - 1; j++)
            {
                int splits = _meshParameters.AbscissaSplits[j];
                double k = _meshParameters.AbscissaK[j] < 0
                    ? -1.0 / _meshParameters.AbscissaK[j]
                    : _meshParameters.AbscissaK[j];
        
                var a = _points[i * totalNx + abscissaSplits];
                var b = _points[i * totalNx + splits + abscissaSplits];
                
                double h = (Math.Abs(k - 1.0) < 1E-14)
                    ? (b.X - a.X) / splits
                    : (b.X - a.X) * (1.0 - k) / (1.0 - Math.Pow(k, splits));
        
                for (int l = 1; l < splits; l++)
                {
                    double x = _points[i * totalNx + abscissaSplits + l - 1].X + h * Math.Pow(k, l - 1);
                    double t = (x - a.X) / (b.X - a.X);
                    double y = a.Y + t * (b.Y - a.Y);
                    
                    _points[i * totalNx + abscissaSplits + l] = new Point(x, y);
                }
        
                abscissaSplits += splits;
            }
        }

        for (int i = 1; i < _meshParameters.AbscissaPointsCount; i++)
        {
            _ix[i] = _ix[i - 1] + _meshParameters.AbscissaSplits[i - 1];
        }
        
        for (int i = 1; i < _meshParameters.OrdinatePointsCount; i++)
        {
            _iy[i] = _iy[i - 1] + _meshParameters.OrdinateSplits[i - 1];
        }

        Points = _points;
    }

    public void CreateElements()
    {
        Span<int> nodes = stackalloc int[4];

        int nx = _meshParameters.AbscissaSplits.Sum() + 1;
        int ny = _meshParameters.OrdinateSplits.Sum() + 1;
        
        for (int i = 0, ielem = 0; i < ny - 1; i++)
        {
            for (int j = 0; j < nx - 1; j++)
            {
                nodes[0] = i * nx + j;
                nodes[1] = i * nx + 1 + j;
                nodes[2] = (i + 1) * nx + j;
                nodes[3] = (i + 1) * nx + 1 + j;
                
                // Find area
                int areaIndex;
                
                for (areaIndex = 0; areaIndex < _meshParameters.Areas.Length; areaIndex++)
                {
                    var area = _meshParameters.Areas[areaIndex];
                    var mxl = _ix[area.LeftBorderNumber];
                    var mxr = _ix[area.RightBorderNumber];
                    var myb = _iy[area.BottomBorderNumber];
                    var myt = _iy[area.TopBorderNumber];

                    if (mxl > j || j + 1 > mxr || myb > i || i + 1 > myt) continue;
                    break;
                }
                
                _elements[ielem++] = new FiniteElement(nodes.ToArray(), _meshParameters.Areas[areaIndex].ParameterNumber);
            }
        }
        
        MarkFictitious();
    }

    private void MarkFictitious()
    {
        _fictitiousElements = Enumerable.Range(0, _elements.Length)
            .Where(ielem => _elements[ielem].AreaNumber == -1)
            .ToList();
        
        foreach (var ielem in _fictitiousElements)
        {
            _fictitiousNodes.AddRange(_elements[ielem].Nodes);
        }

        _fictitiousNodes = _fictitiousNodes
            .GroupBy(node => node)
            .Where(group => group.Count() > 2)
            .SelectMany(group => group)
            .Distinct()
            .ToList();
    }
    
    public void CreateBoundaries()
    {
        foreach (var border in _meshParameters.Borders)
        {
            switch (border.BoundaryType)
            {
                case BoundaryType.Dirichlet:
                    ProcessDirichletBoundary(border);
                    break;
                case BoundaryType.Neumann:
                    ProcessNeumannBoundary(border);
                    break;
                case BoundaryType.Newton:
                    ProcessNewtonBoundary(border);
                    break;
                case BoundaryType.None:
                    break;
                default:
                    ProcessDirichletBoundary(border);
                    break;
            }
        }

        _dirichlet = _dirichlet.DistinctBy(d => d.Node).ToHashSet();
    }

    private void ProcessDirichletBoundary(Border border)
    {
        int nx = _meshParameters.AbscissaPointsCount;
        int totalNx = _meshParameters.AbscissaSplits.Sum() + 1;
        
        for (int i = 0; i < border.PointsIndices.Length - 1; i++)
        {
            int curr = border.PointsIndices[i];
            int next = border.PointsIndices[i + 1];

            int ys = curr / nx;
            int xs = curr - ys * nx;
            int ye = next / nx;
            int xe = next - ye * nx;

            xs = _ix[xs];
            xe = _ix[xe];
            ys = _iy[ys];
            ye = _iy[ye];

            (xs, xe) = xe < xs ? (xe, xs) : (xs, xe);
            (ys, ye) = ye < ys ? (ye, ys) : (ys, ye);

            for (int j = 0; j < _points.Length; j++)
            {
                if (_fictitiousNodes.Any(n => n == j)) continue;

                int iy = j / totalNx;
                int ix = j - iy * totalNx;
                
                if (ix < xs || ix > xe || iy < ys || iy > ye) continue;

                _dirichlet.Add(new Dirichlet(j, _meshParameters.BoundaryFormulas[border.FormulaIndex]));
            }
        }
    }

    private void ProcessNeumannBoundary(Border border)
    {
        _neumann ??= new HashSet<Neumann>();
        
        int nx = _meshParameters.AbscissaPointsCount;
        int totalNx = _meshParameters.AbscissaSplits.Sum() + 1;
        
        for (int i = 0; i < border.PointsIndices.Length - 1; i++)
        {
            int curr = border.PointsIndices[i];
            int next = border.PointsIndices[i + 1];

            int ys = curr / nx;
            int xs = curr - ys * nx;
            int ye = next / nx;
            int xe = next - ye * nx;
            
            xs = _ix[xs];
            xe = _ix[xe];
            ys = _iy[ys];
            ye = _iy[ye];

            (xs, xe) = xe < xs ? (xe, xs) : (xs, xe);
            (ys, ye) = ye < ys ? (ye, ys) : (ys, ye);

            // Horizontal line
            if (ys == ye)
            {
                for (int j = 0; j < _points.Length; j++)
                {
                    int iy = j / totalNx;
                    int ix = j - iy * totalNx;
                    
                    // If point is not on the line
                    if (iy != ys || ix < xs || ix > xe) continue;
                    
                    _neumann.Add(new Neumann(new Edge(j, j + 1),
                        _meshParameters.BoundaryFormulas[border.FormulaIndex]));
                    
                    while (ix + 1 != xe)
                    {
                        j++;
                        ix = j - iy * totalNx;

                        _neumann.Add(new Neumann(new Edge(j, j + 1), 
                            _meshParameters.BoundaryFormulas[border.FormulaIndex]));
                    }
                    break;
                }
            }
            // Vertical line
            else
            {
                for (int j = 0; j < _points.Length; j++)
                {
                    int iy = j / totalNx;
                    int ix = j - iy * totalNx;
                    
                    // If point is not on the line
                    if (ix != xs || iy < ys || iy > ye) continue;
                    
                    _neumann.Add(new Neumann(new Edge(j, j + totalNx),
                        _meshParameters.BoundaryFormulas[border.FormulaIndex]));
                    
                    while (iy + 1 != ye)
                    {
                        j += totalNx;
                        iy = j / totalNx;
                        
                        _neumann.Add(new Neumann(new Edge(j, j + totalNx),
                            _meshParameters.BoundaryFormulas[border.FormulaIndex]));
                    }
                    break;
                }
            }
        }
    }

    private void ProcessNewtonBoundary(Border border)
    {
        _newton ??= [];
        
        int nx = _meshParameters.AbscissaPointsCount;
        int totalNx = _meshParameters.AbscissaSplits.Sum() + 1;
        
        for (int i = 0; i < border.PointsIndices.Length - 1; i++)
        {
            int curr = border.PointsIndices[i];
            int next = border.PointsIndices[i + 1];

            int ys = curr / nx;
            int xs = curr - ys * nx;
            int ye = next / nx;
            int xe = next - ye * nx;
            
            xs = _ix[xs];
            xe = _ix[xe];
            ys = _iy[ys];
            ye = _iy[ye];

            (xs, xe) = xe < xs ? (xe, xs) : (xs, xe);
            (ys, ye) = ye < ys ? (ye, ys) : (ys, ye);

            // Horizontal line
            if (ys == ye)
            {
                for (int j = 0; j < _points.Length; j++)
                {
                    int iy = j / totalNx;
                    int ix = j - iy * totalNx;
                    
                    // If point is not on the line
                    if (iy != ys || ix < xs || ix > xe) continue;

                    _newton.Add(new Newton(new Edge(j, j + 1),
                        _meshParameters.BoundaryFormulas[border.FormulaIndex], 
                        border.FormulaIndex == 2
                            ? _meshParameters.Betas[0]
                            : _meshParameters.Betas[1]));
                    
                    while (ix + 1 != xe)
                    {
                        j++;
                        ix = j - iy * totalNx;

                        _newton.Add(new Newton(new Edge(j, j + 1),
                            _meshParameters.BoundaryFormulas[border.FormulaIndex], 
                            border.FormulaIndex == 2
                                ? _meshParameters.Betas[0]
                                : _meshParameters.Betas[1]));
                    }
                    break;
                }
            }
            // Vertical line
            else
            {
                for (int j = 0; j < _points.Length; j++)
                {
                    int iy = j / totalNx;
                    int ix = j - iy * totalNx;
                    
                    // If point is not on the line
                    if (ix != xs || iy < ys || iy > ye) continue;
                    
                    _newton.Add(new Newton(new Edge(j, j + totalNx),
                        _meshParameters.BoundaryFormulas[border.FormulaIndex], 
                        border.FormulaIndex == 2
                            ? _meshParameters.Betas[0]
                            : _meshParameters.Betas[1]));
                    
                    while (iy + 1 != ye)
                    {
                        j += totalNx;
                        iy = j / totalNx;
                        
                        _newton.Add(new Newton(new Edge(j, j + totalNx),
                            _meshParameters.BoundaryFormulas[border.FormulaIndex], 
                            border.FormulaIndex == 2
                                ? _meshParameters.Betas[0]
                                : _meshParameters.Betas[1]));
                    }
                    break;
                }
            }
        }
    }
    
    public Mesh GetMesh() 
        => new(_points, _elements, _meshParameters.AreaProperties, _dirichlet, _neumann, _newton);
}