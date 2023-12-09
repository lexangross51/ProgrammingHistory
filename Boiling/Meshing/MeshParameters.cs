using Boiling.Meshing.Geometry;
using Boiling.Meshing.JsonConverter;
using Newtonsoft.Json;

namespace Boiling.Meshing;

public struct Area(
    int leftBorderNumber, int rightBorderNumber,
    int bottomBorderNumber, int topBorderNumber,
    int parameterNumber)
{
    public int ParameterNumber { get; set; } = parameterNumber;
    public int LeftBorderNumber { get; set; } = leftBorderNumber;
    public int RightBorderNumber { get; set; } = rightBorderNumber;
    public int BottomBorderNumber { get; set; } = bottomBorderNumber;
    public int TopBorderNumber { get; set; } = topBorderNumber;
}

public class Border(int[] pointsIndices, BoundaryType boundaryType, int formulaIndex)
{
    public int[] PointsIndices { get; set; } = pointsIndices;
    public BoundaryType BoundaryType { get; set; } = boundaryType;
    public int FormulaIndex { get; set; } = formulaIndex;
}

[JsonConverter(typeof(MeshParametersJsonConverter))]
public class MeshParameters
{
    public int AbscissaPointsCount { get; init; }
    public int OrdinatePointsCount { get; init; }
    public Point[] ControlPoints { get; init; } = null!;
    public Area[] Areas { get; init; } = null!;
    public Border[] Borders { get; init; } = null!;
    public Func<double, double, double>[] BoundaryFormulas { get; init; } = null!;
    public double[] Betas { get; init; } = null!;
    public AreaProperty[] AreaProperties { get; init; } = null!;
    public int[] AbscissaSplits { get; init; } = null!;
    public int[] OrdinateSplits { get; init; } = null!;
    public double[] AbscissaK { get; init; } = null!;
    public double[] OrdinateK { get; init; } = null!;
    public int Refinement { get; init; }

    public static MeshParameters ReadJson(string path)
    {
        using var sr = new StreamReader(path);
        return JsonConvert.DeserializeObject<MeshParameters>(sr.ReadToEnd());
    }
}