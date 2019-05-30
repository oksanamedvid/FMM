using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace PFMM
{
    public class MultipoleExpansion
    {
        private static double Distance(Point x, Point y)
        {
            return Math.Sqrt(Math.Pow((x.X - y.X), 2) + Math.Pow((x.Y - y.Y), 2));
        }

        private ulong bcl(int n, int k)
        {
            if (k > n / 2) k = n - k;
            switch (k)
            {
                case 1:
                    return (ulong)n;
                case 0:
                    return 1;
            }

            ulong r;
            if (n + k >= 90)
            {
                // разрядности может не хватить, используем рекурсию
                r = bcl(n - 1, k);
                r += +bcl(n - 1, k - 1);
            }
            else
            {
                r = 1;
                for (int i = 1; i <= k; ++i)
                {
                    r *= (ulong)(n - k + i);
                    r /= (ulong)i;
                }
            }

            return r;
        }

        public Matrix<double> T_ofs(List<Point> points, Point c)
        {
            var tOfs = new double[CustomConstrants.ErrorP, points.Count];
            for (var i = 0; i < CustomConstrants.ErrorP; i++)
            {
                for (var j = 0; j < points.Count; j++)
                {
                    if (i == 0)
                    {
                        tOfs[i, j] = -(1.0);
                    }
                    else
                    {
                        tOfs[i, j] = -(1.0) * -(1.0 / i) * Math.Pow(Distance(points[j], c), i);
                    }
                }
            }

            return DenseMatrix.OfArray(tOfs);
        }

        public Matrix<double> T_ifo(Point c1, Point c2)
        {
            var tIfo = new double[CustomConstrants.ErrorP, CustomConstrants.ErrorP];
            for (var i = 0; i < CustomConstrants.ErrorP; i++)
            {
                if (i == 0)
                {
                    tIfo[i, 0] = Math.Log(Distance(c1, c2));
                }
                else
                {
                    tIfo[i, 0] = -1.0 / (i * Math.Pow(Distance(c1, c2), i));
                }

                for (var j = 1; j < CustomConstrants.ErrorP; j++)
                {
                    if (i == 0)
                    {
                        tIfo[i, j] = Math.Pow(-1, j) * (1.0 / Math.Pow(Distance(c1, c2), j));
                    }
                    else
                    {
                        tIfo[i, j] = Math.Pow(-1, j) * bcl(i + j - 1, j - 1) *
                                     (1.0 / Math.Pow(Distance(c1, c2), i + j));
                    }
                }
            }

            return DenseMatrix.OfArray(tIfo);
        }

        public Matrix<double> T_tfi(List<Point> points, Point c)
        {
            var tIfi = new double[points.Count, CustomConstrants.ErrorP];
            for (var i = 0; i < points.Count; i++)
            {
                for (var j = 0; j < CustomConstrants.ErrorP; j++)
                {
                    tIfi[i, j] = Math.Pow(Distance(points[i], c), j);
                }
            }

            return DenseMatrix.OfArray(tIfi);
        }

        public Matrix<double> T_ofo(Point c1, Point c2)
        {
            var tOfo = new double[CustomConstrants.ErrorP, CustomConstrants.ErrorP];
            for (var i = 0; i < CustomConstrants.ErrorP; i++)
            {
                if (i == 0)
                {
                    tOfo[i, 0] = 1.0;
                }
                else
                {
                    tOfo[i, 0] = -(1.0 / i) * Math.Pow(Distance(c1, c2), i);
                    for (var j = 1; j <= i; j++)
                    {
                        tOfo[i, j] = bcl(i - 1, j - 1) * Math.Pow(Distance(c1, c2), i - j);
                    }
                }
            }

            return DenseMatrix.OfArray(tOfo);
        }

        public Matrix<double> T_ifi(Point c1, Point c2)
        {
            var tIfi = new double[CustomConstrants.ErrorP, CustomConstrants.ErrorP];
            for (var i = 0; i < CustomConstrants.ErrorP; i++)
            {
                for (var j = i; j < CustomConstrants.ErrorP; j++)
                {
                    tIfi[i, j] = bcl(j, i) * Math.Pow(Distance(c1, c2), j - i);
                }
            }

            return DenseMatrix.OfArray(tIfi);
        }
    }
}
