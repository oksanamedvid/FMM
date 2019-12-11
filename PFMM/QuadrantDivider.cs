using System;
using System.Collections.Generic;
namespace PFMM
{
    public static class QuadrantDivider
    {
        public static List<Quadrant> Quadrants = new List<Quadrant>();
        public static double HA = Math.Abs((SquarePoints.A2 - SquarePoints.A1));
        public static double HB = Math.Abs((SquarePoints.B2 - SquarePoints.B1));

        public static void CreateQuadrants(Quadrant parent)
        {
            Quadrants.Clear();
            for (var i = 0; i < 4; ++i)
            {
                Quadrants.Add(new Quadrant(parent));
            }
        }

        public static void InitializeBound(Quadrant quadrant, double a1, double a2, double b1, double b2)
        {
            quadrant.AStart = a1;
            quadrant.AEnd = a2;
            quadrant.BStart = b1;
            quadrant.BEnd = b2;
        }

        public static List<Quadrant> DivideQuadrant(Quadrant quadrant, int l)
        {
            CreateQuadrants(quadrant);
            InitializeBound(Quadrants[0], quadrant.AStart, (HA / Math.Pow(2, l) + quadrant.AStart), quadrant.BStart, (HB / Math.Pow(2, l) + quadrant.BStart));
            InitializeBound(Quadrants[1], (HA / Math.Pow(2, l) + quadrant.AStart), quadrant.AEnd, quadrant.BStart, (HB / Math.Pow(2, l) + quadrant.BStart));
            InitializeBound(Quadrants[2], quadrant.AStart, (HA / Math.Pow(2, l) + quadrant.AStart), (HB / Math.Pow(2, l) + quadrant.BStart), quadrant.BEnd);
            InitializeBound(Quadrants[3], (HA / Math.Pow(2, l) + quadrant.AStart), quadrant.AEnd, (HB / Math.Pow(2, l) + quadrant.BStart), quadrant.BEnd);
            foreach (var point in quadrant.Points)
            {
                if ((point.X >= quadrant.AStart && point.X < (HA / Math.Pow(2, l) + quadrant.AStart)) && (point.Y >= quadrant.BStart && point.Y < (HB / Math.Pow(2, l) + quadrant.BStart)))
                {
                    Quadrants[0].Points.Add(point);
                    continue;
                }

                if ((point.X >= (HA / Math.Pow(2, l) + quadrant.AStart) && point.X <= quadrant.AEnd) && (point.Y >= quadrant.BStart && point.Y < (HB / Math.Pow(2, l) + quadrant.BStart)))
                {
                    Quadrants[1].Points.Add(point);
                    continue;
                }

                if ((point.X >= quadrant.AStart && point.X < (HA / Math.Pow(2, l) + quadrant.AStart)) && (point.Y >= (HB / Math.Pow(2, l) + quadrant.BStart) && point.Y <= quadrant.BEnd))
                {
                    Quadrants[2].Points.Add(point);
                    continue;
                }

                if ((point.X >= (HA / Math.Pow(2, l) + quadrant.AStart) && point.X <= quadrant.AEnd) && (point.Y >= (HB / Math.Pow(2, l) + quadrant.BStart) && point.Y <= quadrant.BEnd))
                {
                    Quadrants[3].Points.Add(point);
                }
            }

            return Quadrants;
        }

        public static List<Quadrant> NearField(Quadrant quad, List<Quadrant> quadrantsOfLevel, int l)
        {
            var nearField = new List<Quadrant>();
            var tempA = HA / Math.Pow(2, l);
            var tempB = HB / Math.Pow(2, l);

            foreach (var quadrant in quadrantsOfLevel)
            {
                if ((Math.Abs(Math.Abs(quad.AStart - quadrant.AStart) - (tempA)) < Double.Epsilon ||
                     (Math.Abs(Math.Abs(quad.AStart - quadrant.AStart)) < Double.Epsilon)) &&
                    ((Math.Abs(Math.Abs(quad.BStart - quadrant.BStart) - (tempB)) < Double.Epsilon) ||
                     (Math.Abs(Math.Abs(quad.BStart - quadrant.BStart)) < Double.Epsilon)))
                {
                    nearField.Add(quadrant);
                }
            }

            return nearField;
        }

        public static List<Quadrant> FarField(Quadrant quad, List<Quadrant> quadrantsOfLevel, int l)
        {
            var farField = new List<Quadrant>();
            var tempA = HA / Math.Pow(2, l);
            var tempB = HB / Math.Pow(2, l);

            foreach (var quadrant in quadrantsOfLevel)
            {
                if (Math.Abs(quad.AStart - quadrant.AStart) > tempA || Math.Abs(quad.BStart - quadrant.BStart) > tempB)
                {
                    farField.Add(quadrant);
                }
            }

            return farField;
        }
    }
}
