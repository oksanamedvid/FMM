using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace PFMM
{
    public class FastMultipoleMethod
    {
        public int Level;

        public Dictionary<Point, double> MethodsMultLevel(List<Point> points)
        {
            Level = 1;
            var wholeQuadrant = new Quadrant(points);
            var currentLevelOfQuadrants = new List<Quadrant>();

            var quadrantsOfLevels = new List<List<Quadrant>>
            {
                new List<Quadrant>()
            };

            AssignListQuadrants(currentLevelOfQuadrants, QuadrantDivider.DivideQuadrant(wholeQuadrant, Level));
            AssignListQuadrants(quadrantsOfLevels[Level - 1], currentLevelOfQuadrants);
            while (MaxPointCountInQuadrant(currentLevelOfQuadrants) > CustomConstants.MaxPointCountInQuad)
            {
                Level++;
                quadrantsOfLevels.Add(new List<Quadrant>());
                AssignListQuadrants(currentLevelOfQuadrants, NextLevel(currentLevelOfQuadrants));
                AssignListQuadrants(quadrantsOfLevels[Level - 1], currentLevelOfQuadrants);
            }

            MultipleLevelsInteraction(quadrantsOfLevels);

            var dictionary = new Dictionary<Point, double>();
            foreach (var quad in quadrantsOfLevels[Level - 1])
            {
                for (var j = 0; j < quad.Points.Count; j++)
                {
                    dictionary.Add(quad.Points[j], quad.Interaction[j]);
                }
            }

            return dictionary;
        }

        public Dictionary<Point, double> MethodsSingleLevel(List<Point> points)
        {
            Level = 1;
            var wholeQuadrant = new Quadrant(points);
            var currentLevelOfQuadrants = new List<Quadrant>();
            var quadrantsOfLevels = new List<List<Quadrant>>
            {
                new List<Quadrant>()
            };

            AssignListQuadrants(currentLevelOfQuadrants, QuadrantDivider.DivideQuadrant(wholeQuadrant, Level));
            AssignListQuadrants(quadrantsOfLevels[Level - 1], currentLevelOfQuadrants);

            while (MaxPointCountInQuadrant(currentLevelOfQuadrants) > CustomConstants.MaxPointCountInQuad)
            {
                Level++;
                quadrantsOfLevels.Add(new List<Quadrant>());
                AssignListQuadrants(currentLevelOfQuadrants, NextLevel(currentLevelOfQuadrants));
                AssignListQuadrants(quadrantsOfLevels[Level - 1], currentLevelOfQuadrants);
            }

            SingleLevelInteraction(quadrantsOfLevels);

            var dictionary = new Dictionary<Point, double>();
            foreach (var quad in quadrantsOfLevels[Level - 1])
            {
                for (var j = 0; j < quad.Points.Count; j++)
                {
                    dictionary.Add(quad.Points[j], quad.Interaction[j]);
                }
            }

            return dictionary;
        }

        private void MultipleLevelsInteraction(IReadOnlyList<List<Quadrant>> quadrantsOfLevels)
        {
            var multEx = new MultipoleExpansion();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            for (var i = 0; i < quadrantsOfLevels[Level - 1].Count; i++)
            {
                if (quadrantsOfLevels[Level - 1][i]?.Points.Any() != true)
                {
                    continue;
                }

                var weight = new List<double>();
                quadrantsOfLevels[Level - 1][i].Points.ForEach(p => { weight.Add(p.AdditionalValue); });
                quadrantsOfLevels[Level - 1][i].OutgoingExpansions = MultiplyMatrixOnVector(multEx.T_ofs(quadrantsOfLevels[Level - 1][i].Points,
                        quadrantsOfLevels[Level - 1][i].CentreQuadrant()), DenseVector.OfArray(weight.ToArray()));
            }
            Console.WriteLine(watch.ElapsedMilliseconds);
            watch.Stop();
            watch = System.Diagnostics.Stopwatch.StartNew();

            var child = new List<Quadrant>();
            for (var i = Level - 2; i > -1; --i)
            {
                for (var j = 0; j < quadrantsOfLevels[i].Count; ++j)
                {
                    quadrantsOfLevels[i][j].OutgoingExpansions = DenseVector.OfArray(new double[CustomConstants.ErrorP]);
                    if (quadrantsOfLevels[i][j]?.Points.Any() != true)
                    {
                        continue;
                    }

                    AssignListQuadrants(child, GetChildren(quadrantsOfLevels[i][j], quadrantsOfLevels[i + 1]));
                    foreach (var quad in child)
                    {
                        if (quad.Points?.Any() != true)
                        {
                            continue;
                        }

                        quadrantsOfLevels[i][j].OutgoingExpansions += MultiplyMatrixOnVector(
                            multEx.T_ofo(quadrantsOfLevels[i][j].CentreQuadrant(), quad.CentreQuadrant()),
                            quad.OutgoingExpansions);
                    }
                }
            }

            Console.WriteLine(watch.ElapsedMilliseconds);
            watch.Stop();
            watch = System.Diagnostics.Stopwatch.StartNew();

            var interactionList = new List<Quadrant>();
            var temp = new List<Quadrant>();
            var parentNearField = new List<Quadrant>();
            var childrenOfParentNeighbors = new List<Quadrant>();

            for (var j = 0; j < quadrantsOfLevels[0].Count; ++j)
            {
                quadrantsOfLevels[0][j].IncomingExpansions = DenseVector.OfArray(new double[CustomConstants.ErrorP]);
            }

            for (var i = 1; i < Level; i++)
            {
                for (var j = 0; j < quadrantsOfLevels[i].Count; ++j)
                {
                    quadrantsOfLevels[i][j].IncomingExpansions = DenseVector.OfArray(new double[CustomConstants.ErrorP]);

                    if (quadrantsOfLevels[i][j].Points?.Any() != true)
                    {
                        continue;
                    }

                    AssignListQuadrants(temp, quadrantsOfLevels[i - 1]);
                    AssignListQuadrants(parentNearField, QuadrantDivider.NearField(quadrantsOfLevels[i][j].Parent, temp, i));

                    var theSameQuad = parentNearField.FirstOrDefault(q =>
                        Math.Abs(quadrantsOfLevels[i][j].Parent.AEnd - q.AEnd) < Double.Epsilon
                        && Math.Abs(quadrantsOfLevels[i][j].Parent.BEnd - q.BEnd) < Double.Epsilon
                        && Math.Abs(quadrantsOfLevels[i][j].Parent.BStart - q.BStart) < Double.Epsilon
                        && Math.Abs(quadrantsOfLevels[i][j].Parent.AStart - q.AStart) < Double.Epsilon);
                    parentNearField.Remove(theSameQuad);

                    childrenOfParentNeighbors.Clear();
                    foreach (var nearFild in parentNearField)
                    {
                        childrenOfParentNeighbors.AddRange(GetChildren(nearFild, quadrantsOfLevels[i]));
                    }

                    AssignListQuadrants(interactionList,
                        QuadrantDivider.FarField(quadrantsOfLevels[i][j], childrenOfParentNeighbors, i + 1));

                    foreach (var quad in interactionList)
                    {
                        if (quad.Points == null || !quad.Points.Any())
                        {
                            continue;
                        }

                        quadrantsOfLevels[i][j].IncomingExpansions += MultiplyMatrixOnVector(multEx.T_ifo(
                            quadrantsOfLevels[i][j].CentreQuadrant(),
                            quad.CentreQuadrant()), quad.OutgoingExpansions);
                    }
                }
            }

            Console.WriteLine(watch.ElapsedMilliseconds);
            watch.Stop();
            watch = System.Diagnostics.Stopwatch.StartNew();

            for (var i = 1; i < Level; ++i)
            {
                for (var j = 0; j < quadrantsOfLevels[i].Count; ++j)
                {
                    if (quadrantsOfLevels[i][j].Points?.Any() != true)
                    {
                        continue;
                    }

                    var parent = GetParentQuad(quadrantsOfLevels[i][j].Parent, quadrantsOfLevels[i - 1]);
                    var quadCentre = quadrantsOfLevels[i][j].CentreQuadrant();
                    var parentCentre = parent.CentreQuadrant();
                    quadrantsOfLevels[i][j].IncomingExpansions +=
                        MultiplyMatrixOnVector(multEx.T_ifi(quadCentre, parentCentre), parent.IncomingExpansions);
                }
            }

            Console.WriteLine(watch.ElapsedMilliseconds);
            watch.Stop();
            watch = System.Diagnostics.Stopwatch.StartNew();

            var nearList = new List<Quadrant>();
            for (var j = 0; j < quadrantsOfLevels[Level - 1].Count; ++j)
            {
                if (quadrantsOfLevels[Level - 1][j].Points == null || !quadrantsOfLevels[Level - 1][j].Points.Any())
                {
                    continue;
                }

                var a = MultiplyMatrixOnVector(
                    multEx.T_tfi(quadrantsOfLevels[Level - 1][j].Points,
                        quadrantsOfLevels[Level - 1][j].CentreQuadrant()),
                    quadrantsOfLevels[Level - 1][j].IncomingExpansions);

                quadrantsOfLevels[Level - 1][j].Interaction = a;

                AssignListQuadrants(nearList, QuadrantDivider.NearField(quadrantsOfLevels[Level - 1][j], temp, Level));

                foreach (var quad in nearList)
                {
                    if (quad.Points == null || !quad.Points.Any())
                    {
                        continue;
                    }

                    List<double> weightNearQuad = new List<double>();
                    quad.Points.ForEach(p => { weightNearQuad.Add(p.AdditionalValue); });

                    var c = MultiplyMatrixOnVector(A(quadrantsOfLevels[Level - 1][j].Points, quad.Points),
                        DenseVector.OfArray(weightNearQuad.ToArray()));

                    quadrantsOfLevels[Level - 1][j].Interaction += c;
                }
            }

            Console.WriteLine(watch.ElapsedMilliseconds);
            watch.Stop();
        }

        public Dictionary<Point, double> SolveDirihleProblem(List<Point> pointsToSolve, List<Point> pointOnBounds)
        {
            var solveQuad = new Quadrant(pointsToSolve);
            var boundsQuad = new Quadrant(pointOnBounds);
            boundsQuad.AStart = -2 ;
            boundsQuad.AEnd = 6;
            boundsQuad.BStart = -4 ;
            boundsQuad.BEnd = 4;

            var multEx = new MultipoleExpansion();

            var weight = new List<double>();
            boundsQuad.Points.ForEach(p => { weight.Add(p.AdditionalValue); });

            boundsQuad.OutgoingExpansions = MultiplyMatrixOnVector(multEx.T_ofs(boundsQuad.Points,
                boundsQuad.CentreQuadrant()), DenseVector.OfArray(weight.ToArray()));

            solveQuad.IncomingExpansions = MultiplyMatrixOnVector(multEx.T_ifo(
                solveQuad.CentreQuadrant(), boundsQuad.CentreQuadrant()), boundsQuad.OutgoingExpansions);

            solveQuad.Interaction = MultiplyMatrixOnVector(multEx.T_tfi(solveQuad.Points,
                    solveQuad.CentreQuadrant()),
                solveQuad.IncomingExpansions);

            var dictionary = new Dictionary<Point, double>();
            for (var j = 0; j < solveQuad.Points.Count; j++)
            {
                dictionary.Add(solveQuad.Points[j], solveQuad.Interaction[j]);
            }

            return dictionary;
        }

        private void SingleLevelInteraction(IReadOnlyList<List<Quadrant>> quadrantsOfLevels)
        {
            var multEx = new MultipoleExpansion();
            for (var i = 0; i < quadrantsOfLevels[Level - 1].Count; i++)
            {
                if (quadrantsOfLevels[Level - 1][i].Points?.Any() != true)
                {
                    continue;
                }

                var weight = new List<double>();
                quadrantsOfLevels[Level - 1][i].Points.ForEach(p => { weight.Add(p.AdditionalValue); });

                quadrantsOfLevels[Level - 1][i].OutgoingExpansions
                    = MultiplyMatrixOnVector(multEx.T_ofs(quadrantsOfLevels[Level - 1][i].Points,
                        quadrantsOfLevels[Level - 1][i].CentreQuadrant()), DenseVector.OfArray(weight.ToArray()));
            }

            var interactionList = new List<Quadrant>();
            var temp = new List<Quadrant>();
            for (var i = 0; i < quadrantsOfLevels[Level - 1].Count; ++i)
            {
                if (quadrantsOfLevels[Level - 1][i].Points?.Any() != true)
                {
                    continue;
                }

                AssignListQuadrants(temp, quadrantsOfLevels[Level - 1]);
                AssignListQuadrants(interactionList, QuadrantDivider.FarField(quadrantsOfLevels[Level - 1][i], temp, Level));

                quadrantsOfLevels[Level - 1][i].IncomingExpansions = DenseVector.OfArray(new double[CustomConstants.ErrorP]);
                foreach (var quad in interactionList)
                {
                    if (quad.Points?.Any() != true)
                    {
                        continue;
                    }

                    quadrantsOfLevels[Level - 1][i].IncomingExpansions += MultiplyMatrixOnVector(multEx.T_ifo(
                        quadrantsOfLevels[Level - 1][i].CentreQuadrant(), quad.CentreQuadrant()), quad.OutgoingExpansions);
                }
            }

            var nearList = new List<Quadrant>();
            for (var j = 0; j < quadrantsOfLevels[Level - 1].Count; ++j)
            {
                if (quadrantsOfLevels[Level - 1][j].Points?.Any() != true)
                {
                    continue;
                }

                var a = MultiplyMatrixOnVector( multEx.T_tfi(quadrantsOfLevels[Level - 1][j].Points,
                        quadrantsOfLevels[Level - 1][j].CentreQuadrant()),
                    quadrantsOfLevels[Level - 1][j].IncomingExpansions);

                quadrantsOfLevels[Level - 1][j].Interaction = a;

                AssignListQuadrants(nearList, QuadrantDivider.NearField(quadrantsOfLevels[Level - 1][j], temp, Level));

                foreach (var quad in nearList)
                {
                    if (quad.Points?.Any() != true)
                    {
                        continue;
                    }

                    List<double> weightNearQuad = new List<double>();
                    quad.Points.ForEach(p => { weightNearQuad.Add(p.AdditionalValue); });

                    var c = MultiplyMatrixOnVector(A(quadrantsOfLevels[Level - 1][j].Points, quad.Points),
                        DenseVector.OfArray(weightNearQuad.ToArray()));

                    quadrantsOfLevels[Level - 1][j].Interaction += c;
                }
            }
        }

        private Vector<double> MultiplyMatrixOnVector(Matrix<double> matrix, Vector<double> vector)
        {
            var resultVector = new double[matrix.RowCount];

            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    resultVector[i] += matrix[i, j] * vector[j];
                }
            }
            return DenseVector.OfArray(resultVector);
        }

        private Quadrant GetParentQuad(Quadrant par, IEnumerable<Quadrant> quadrants)
        {
            return quadrants.FirstOrDefault(q =>
                Math.Abs(q.AStart - par.AStart) < Double.Epsilon && Math.Abs(q.AEnd - par.AEnd) < Double.Epsilon &&
                Math.Abs(q.BStart - par.BStart) < Double.Epsilon && Math.Abs(q.BEnd - par.BEnd) < Double.Epsilon);
        }

        private Matrix<double> A(IReadOnlyList<Point> tauPoints, IReadOnlyList<Point> tPoints)
        {
            var a = new double[tauPoints.Count, tPoints.Count];
            for (var i = 0; i < tauPoints.Count; i++)
            {
                for (var j = 0; j < tPoints.Count; j++)
                {

                    a[i, j] = Math.Abs(tauPoints[i].X - tPoints[j].X) > Double.Epsilon &&
                              Math.Abs(tauPoints[i].Y - tPoints[j].Y) > Double.Epsilon
                        ? NearInteraction(tauPoints[i], tPoints[j])
                        : 0;
                }
            }

            return DenseMatrix.OfArray(a);
        }

        private double NearInteraction(Point p1, Point p2)
        {
            return Math.Log(Math.Sqrt(Math.Pow((p1.X - p2.X), 2) + Math.Pow((p1.Y - p2.Y), 2)));
        }

        public void InitializeListQuadrants(List<Quadrant> partOfQuadrants, List<Quadrant> quadrants)
        {
            quadrants.InsertRange(quadrants.Count, partOfQuadrants);
        }

        public void AssignQuadrant(Quadrant assignTo, Quadrant assignFrom)
        {
            AssignPoint(assignTo.Points, assignFrom.Points);
            assignTo.AEnd = assignFrom.AEnd;
            assignTo.BEnd = assignFrom.BEnd;
            assignTo.AStart = assignFrom.AStart;
            assignTo.BStart = assignFrom.BStart;
        }

        private void AssignListQuadrants(List<Quadrant> assignTo, List<Quadrant> assignFrom)
        {
            assignTo.Clear();
            assignTo.InsertRange(0, assignFrom);
        }

        private void AssignPoint(List<Point> assignTo, List<Point> assignFrom)
        {
            assignTo.Clear();
            assignTo.InsertRange(0, assignFrom);
        }

        private int MaxPointCountInQuadrant(List<Quadrant> quadrantOfCurrentLevel)
        {
            var max = 0;
            foreach (var quad in quadrantOfCurrentLevel)
            {
                if (quad.Points.Count > max)
                {
                    max = quad.Points.Count;
                }
            }

            return max;
        }

        private List<Quadrant> NextLevel(IEnumerable<Quadrant> quadrantsOfCurrentLevel)
        {
            var childrenOfQuadrant = new List<Quadrant>();
            var childrenOfAllQuadrants = new List<Quadrant>();
            var currentQuadrant = new Quadrant();
            foreach (var quad in quadrantsOfCurrentLevel)
            {
                AssignQuadrant(currentQuadrant, quad);
                if (!currentQuadrant.Points.Any())
                {
                    continue;
                }
                AssignListQuadrants(childrenOfQuadrant, QuadrantDivider.DivideQuadrant(currentQuadrant, Level));
                InitializeListQuadrants(childrenOfQuadrant, childrenOfAllQuadrants);
            }

            return childrenOfAllQuadrants;
        }

        private List<Quadrant> GetChildren(Quadrant parent, IEnumerable<Quadrant> quadrants)
        {
            var copyOfQuadrants = new List<Quadrant>(quadrants);
            return copyOfQuadrants.Where(quadrant => Math.Abs(quadrant.Parent.AEnd - parent.AEnd) < Double.Epsilon
                                               && Math.Abs(quadrant.Parent.BEnd - parent.BEnd) < Double.Epsilon
                                               && Math.Abs(quadrant.Parent.AStart - parent.AStart) < Double.Epsilon
                                               && Math.Abs(quadrant.Parent.BStart - parent.BStart) < Double.Epsilon).ToList();
        }
    }
}