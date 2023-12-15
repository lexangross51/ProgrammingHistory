using System.Globalization;
using System.Linq.Expressions;
using Boiling.Meshing;
using Boiling.Meshing.Geometry;

namespace Boiling;

public static class Utilities
{
    public static Func<double, double, double> GetFunctionFromString(string expression)
    {
        var x = Expression.Parameter(typeof(double), "x");
        var y = Expression.Parameter(typeof(double), "y");
        var body = System.Linq.Dynamic.Core.DynamicExpressionParser.ParseLambda(new[] { x, y }, null, expression).Body;

        return Expression.Lambda<Func<double, double, double>>(body, x, y).Compile();
    }

    public static IEnumerable<Point> GeneratePoints(double xMin, double xMax, double yMin, double yMax, int nx, int ny)
    {
        var points = new List<Point>(nx * ny);

        double hx = (xMax - xMin) / (nx - 1);
        double hy = (yMax - yMin) / (ny - 1);
        
        for (int i = 0; i < ny; i++)
        {
            for (int j = 0; j < nx; j++)
            {
                points.Add(new Point(xMin + j * hx, yMin + i * hy));
            }
        }
        
        return points;
    }
    
    public static void SaveMesh(Mesh mesh, string folder)
    {
        // Points
        var sw = new StreamWriter($"{folder}/points");
        foreach (var p in mesh.Points)
        {
            sw.WriteLine($"{p.X} {p.Y}", CultureInfo.CurrentCulture);
        }
        sw.Close();

        // Elements
        sw = new StreamWriter($"{folder}/elements");
        foreach (var element in mesh.Elements)
        {
            var nodes = element.Nodes;
            sw.WriteLine($"{nodes[0]} {nodes[1]} {nodes[3]} {nodes[2]}");
        }
        sw.Close();
        
        // Dirichlet
        sw = new StreamWriter($"{folder}/dirichlet");
        foreach (var dir in mesh.Dirichlet)
        {
            sw.WriteLine($"{dir.Node} {dir.Value}");
        }
        sw.Close();
        
        if (mesh.Neumann is null) return;
        
        // Neumann
        sw = new StreamWriter($"{folder}/neumann");
        foreach (var n in mesh.Neumann)
        {
            var p1 = mesh.Points[n.Border.Node1];
            var p2 = mesh.Points[n.Border.Node2];
            
            sw.WriteLine($"{n.Border.Node1} {n.Border.Node2} {n.Theta(p1.X, p1.Y)} {n.Theta(p2.X, p2.Y)}");
        }
        sw.Close();
    }

    public static void ElementsPerAreas(Mesh mesh, string folder)
    {
        var p = mesh.Points;
        
        using var sw = new StreamWriter($"{folder}/elemsAreas");

        foreach (var e in mesh.Elements)
        {
            var n = e.Nodes;
            sw.WriteLine($"{p[n[0]].X} {p[n[0]].Y} {p[n[^1]].X} {p[n[^1]].Y} {e.VelocityArea}");
        }
    }

    public static void PrintAtTime(Mesh mesh, double[] solution, int timeMoment, string folder)
    {
        using var sw = new StreamWriter($"{folder}/solution{timeMoment}");
        
        for (int i = 0; i < solution.Length; i++)
        {
            var p = mesh.Points[i];
            sw.WriteLine($"{p.X} {p.Y} {solution[i]}");
        }
    }

    public static void PrintVelocities(Mesh mesh, Func<int, Point> vFunc, string folder)
    {
        var wx = new StreamWriter($"{folder}/velocityX");
        for (int ielem = 0; ielem < mesh.Elements.Length; ielem++)
        {
            var nodes = mesh.Elements[ielem].Nodes;
            double x = nodes.Select(i => mesh.Points[i].X).Sum() / 4.0;
            double y = nodes.Select(i => mesh.Points[i].Y).Sum() / 4.0;
            var v = vFunc(ielem);
            
            wx.WriteLine($"{x} {y} {v.X}");
        }
        wx.Close();
        
        var wy = new StreamWriter($"{folder}/velocityY");
        for (int ielem = 0; ielem < mesh.Elements.Length; ielem++)
        {
            var nodes = mesh.Elements[ielem].Nodes;
            double x = nodes.Select(i => mesh.Points[i].X).Sum() / 4.0;
            double y = nodes.Select(i => mesh.Points[i].Y).Sum() / 4.0;
            var v = vFunc(ielem);
            
            wy.WriteLine($"{x} {y} {v.Y}");
        }
        wy.Close();
    }
}

public static class EnumerableExtensions
{
    public static double Norm(this IEnumerable<double> collection)
        => Math.Sqrt(collection.Sum(item => item * item));

    public static double Dot(this double[] first, double[] second)
        => first.Select((t, i) => t * second[i]).Sum();
}