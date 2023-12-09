namespace Boiling.Meshing.Geometry;

public struct Interval(double left, double right)
{
    public double LeftBorder { get; set; } = left;
    public double RightBorder { get; set; } = right;

    public double Length => RightBorder - LeftBorder;
}

public struct Rectangle(Point leftBottom, Point rightTop)
{
    public Point LeftBottom { get; set; } = leftBottom;
    public Point RightTop { get; set; } = rightTop;

    public static double Square(Point leftBottom, Point rightTop) 
        => (rightTop.X - leftBottom.X) * (rightTop.Y - leftBottom.Y);
}