using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;

namespace PFMM
{
    public class Quadrant
    {
        public List<Point> Points;
        public double AStart;
        public double AEnd;
        public double BStart;
        public double BEnd;
        public Quadrant Parent;

        public Vector<double> OutgoingExpansions;
        public Vector<double> IncomingExpansions;
        public Vector<double> Interaction;

        public Quadrant(List<Point> p)
        {
            Points = new List<Point>();
            Parent = new Quadrant();
            AssignPoints(p);
            AStart = SquarePoints.A1;
            AEnd = SquarePoints.A2;
            BStart = SquarePoints.B1;
            BEnd = SquarePoints.B2;
        }

        public Quadrant()
        {
            Points = new List<Point>();
            AStart = SquarePoints.A1;
            AEnd = SquarePoints.A2;
            BStart = SquarePoints.B1;
            BEnd = SquarePoints.B2;
        }

        public Quadrant(Quadrant p)
        {
            Points = new List<Point>();
            Parent = new Quadrant();
            AssignQuadrant(p);
            AStart = SquarePoints.A1;
            AEnd = SquarePoints.A2;
            BStart = SquarePoints.B1;
            BEnd = SquarePoints.B2;
        }

        public void AssignPoints(List<Point> p)
        {
            foreach (var point in p)
            {
                Points.Add(point);
            }
        }

        public void AssignQuadrant(Quadrant p)
        {
            Parent.AEnd = p.AEnd;
            Parent.BEnd = p.BEnd;
            Parent.AStart = p.AStart;
            Parent.BStart = p.BStart;
            Parent.Points = new List<Point>();
            for (var i = 0; i < p.Points.Count; ++i)
            {
                Parent.Points.Add(p.Points[i]);
            }
        }

        public Point CentreQuadrant()
        {
            return new Point
            {
                X = AStart + (AEnd - AStart) / 2.0,
                Y = BStart + (BEnd - BStart) / 2.0
            };
        }
    }
}
