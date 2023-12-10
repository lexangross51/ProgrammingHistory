namespace Boiling.Meshing;

public struct AreaProperty(double lambda, double rhoCp, 
    Func<double, double, double> source)
{
    public double Lambda { get; set; } = lambda;
    public double RhoCp { get; set; } = rhoCp;
    public Func<double, double, double> Source { get; set; } = source;
}