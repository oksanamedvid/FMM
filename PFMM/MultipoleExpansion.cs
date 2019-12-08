using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace PFMM
{
    public class MultipoleExpansion
    {
        private double Distance(Point x, Point y)
        {
            return Math.Sqrt(Math.Pow((x.X - y.X), 2) + Math.Pow((x.Y - y.Y), 2));
        }

        public double Factorial(int number)
        {
            if (number < 2)
                return 1;
            return number * Factorial(number - 1);
        }

        private double Bcl(int n, int k)
        {
            return (Factorial(n)) / (Factorial(k) * Factorial(n - k));
        }

        public Matrix<double> T_ofs(List<Point> points, Point c)
        {
            var tOfs = new double[CustomConstants.ErrorP, points.Count];
            for (var i = 0; i < CustomConstants.ErrorP; i++)
            {
                for (var j = 0; j < points.Count; j++)
                {
                    if (i == 0)
                    {
                        tOfs[i, j] = (1.0);
                    }
                    else
                    {
                        tOfs[i, j] = (-1.0 / i) * Math.Pow(Distance(points[j], c), i);
                    }
                }
            }

            return DenseMatrix.OfArray(tOfs);
        }

        public Matrix<double> T_ifo(Point c1, Point c2)
        {
            var tIfo = new double[CustomConstants.ErrorP, CustomConstants.ErrorP];
            for (var i = 0; i < CustomConstants.ErrorP; i++)
            {
                if (i == 0)
                {
                    tIfo[i, 0] = Math.Log(Distance(c1, c2));
                }
                else
                {
                    tIfo[i, 0] = (-1.0 / (i * Math.Pow(Distance(c1, c2), i)));
                }

                for (var j = 1; j < CustomConstants.ErrorP; j++)
                {
                    if (i == 0)
                    {
                        tIfo[i, j] = Math.Pow(-1, j) * (1.0 / Math.Pow(Distance(c1, c2), j));
                    }
                    else
                    {
                        tIfo[i, j] = Math.Pow(-1, j) * Bcl(i + j - 1, j - 1) *
                                     (1.0 / Math.Pow(Distance(c1, c2), i + j));
                    }
                }
            }

            return DenseMatrix.OfArray(tIfo);
        }

        public Matrix<double> T_tfi(List<Point> points, Point c)
        {
            var tIfi = new double[points.Count, CustomConstants.ErrorP];
            for (var i = 0; i < points.Count; i++)
            {
                for (var j = 0; j < CustomConstants.ErrorP; j++)
                {
                    tIfi[i, j] = Math.Pow(Distance(points[i], c), j);
                }
            }

            return DenseMatrix.OfArray(tIfi);
        }

        public Matrix<double> T_ofo(Point c1, Point c2)
        {
            var tOfo = new double[CustomConstants.ErrorP, CustomConstants.ErrorP];
            for (var i = 0; i < CustomConstants.ErrorP; i++)
            {
                if (i == 0)
                {
                    tOfo[0, 0] = 1.0;
                }
                else
                {
                    tOfo[i, 0] = -(1.0 / i) * Math.Pow(Distance(c1, c2), i);
                    for (var j = 1; j <= i; j++)
                    {
                        tOfo[i, j] = Bcl(i - 1, j - 1) * Math.Pow(Distance(c1, c2), i - j);
                    }
                }
            }

            return DenseMatrix.OfArray(tOfo);
        }

        public Matrix<double> T_ifi(Point c1, Point c2)
        {
            var tIfi = new double[CustomConstants.ErrorP, CustomConstants.ErrorP];
            for (var i = 0; i < CustomConstants.ErrorP; i++)
            {
                for (var j = i; j < CustomConstants.ErrorP; j++)
                {
                    tIfi[i, j] = Bcl(j, i) * Math.Pow(Distance(c1, c2), j - i);
                }
            }

            return DenseMatrix.OfArray(tIfi);
        }
    }
}
