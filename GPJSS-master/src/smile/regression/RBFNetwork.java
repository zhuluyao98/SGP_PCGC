/*******************************************************************************
 * Copyright (c) 2010 Haifeng Li
 *   
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *  
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *******************************************************************************/

package smile.regression;

import smile.math.distance.EuclideanDistance;
import smile.math.distance.Metric;
import smile.math.matrix.DenseMatrix;
import smile.math.matrix.Matrix;
import smile.math.matrix.QR;
import smile.math.rbf.GaussianRadialBasis;
import smile.math.rbf.RadialBasisFunction;
import smile.util.SmileUtils;

import java.util.Arrays;

/**
 * Radial basis function network. A radial basis function network is an
 * artificial neural network that uses radial basis functions as activation
 * functions. It is a linear combination of radial basis functions. They are
 * used in function approximation, time series prediction, and control.
 * <p>
 * In its basic form, radial basis function network is in the form
 * <p>
 * y(x) = &Sigma; w<sub>i</sub> &phi;(||x-c<sub>i</sub>||)
 * <p>
 * where the approximating function y(x) is represented as a sum of N radial
 * basis functions &phi;, each associated with a different center c<sub>i</sub>,
 * and weighted by an appropriate coefficient w<sub>i</sub>. For distance,
 * one usually chooses Euclidean distance. The weights w<sub>i</sub> can
 * be estimated using the matrix methods of linear least squares, because
 * the approximating function is linear in the weights.
 * <p>
 * The points c<sub>i</sub> are often called the centers of the RBF networks,
 * which can be randomly selected from training data, or learned by some clustering
 * method (e.g. k-means), or learned together with weight parameters undergo
 * a supervised learning processing (e.g. error-correction learning).
 * <p>
 * Popular choices for &phi; comprise the Gaussian function and the so
 * called thin plate splines. The advantage of the thin plate splines is that
 * their conditioning is invariant under scalings. Gaussian, multi-quadric
 * and inverse multi-quadric are infinitely smooth and and involve a scale
 * or shape parameter, r<sub><small>0</small></sub> &gt; 0. Decreasing
 * r<sub><small>0</small></sub> tends to flatten the basis function. For a
 * given function, the quality of approximation may strongly depend on this
 * parameter. In particular, increasing r<sub><small>0</small></sub> has the
 * effect of better conditioning (the separation distance of the scaled points
 * increases).
 * <p>
 * A variant on RBF networks is normalized radial basis function (NRBF)
 * networks, in which we require the sum of the basis functions to be unity.
 * NRBF arises more naturally from a Bayesian statistical perspective. However,
 * there is no evidence that either the NRBF method is consistently superior
 * to the RBF method, or vice versa.
 *
 * <h2>References</h2>
 * <ol>
 * <li> Simon Haykin. Neural Networks: A Comprehensive Foundation (2nd edition). 1999. </li> 
 * <li> T. Poggio and F. Girosi. Networks for approximation and learning. Proc. IEEE 78(9):1484-1487, 1990. </li>
 * <li> Nabil Benoudjit and Michel Verleysen. On the kernel widths in radial-basis function networks. Neural Process, 2003.</li>
 * </ol>
 * 
 * @see RadialBasisFunction
 * @see SVR
 * 
 * @author Haifeng Li
 */
public class RBFNetwork<T> implements Regression<T> {
    private static final long serialVersionUID = 1L;

    /**
     * The centers of RBF functions.
     */
    private T[] centers;
    /**
     * The linear weights.
     */
    private double[] w;
    /**
     * The distance functor.
     */
    private Metric<T> distance;
    /**
     * The radial basis functions.
     */
    private RadialBasisFunction[] rbf;
    /**
     * True to fit a normalized RBF network.
     */
    private boolean normalized;

    /**
     * Trainer for RBF networks.
     */
    public static class Trainer<T> extends RegressionTrainer<T> {
        /**
         * The number of centers.
         */
        private int m = 10;
        /**
         * The distance metric functor.
         */
        private Metric<T> distance;
        /**
         * The radial basis functions.
         */
        private RadialBasisFunction[] rbf;
        /**
         * True to fit a normalized RBF network.
         */
        private boolean normalized = false;

        /**
         * Constructor.
         * 
         * @param distance the distance metric functor.
         */
        public Trainer(Metric<T> distance) {
            this.distance = distance;
        }
        
        /**
         * Sets the radial basis function.
         * @param rbf the radial basis function.
         * @param m the number of basis functions.
         */
        public Trainer<T> setRBF(RadialBasisFunction rbf, int m) {
            this.m = m;
            this.rbf = rep(rbf, m);
            return this;
        }
        
        /**
         * Sets the radial basis functions.
         * @param rbf the radial basis functions.
         */
        public Trainer<T> setRBF(RadialBasisFunction[] rbf) {
            this.m = rbf.length;
            this.rbf = rbf;
            return this;
        }
        
        /**
         * Sets the number of centers.
         * @param m the number of centers.
         */
        public Trainer<T> setNumCenters(int m) {
            this.m = m;
            return this;
        }
        
        /**
         * Sets true to learn normalized RBF network.
         * @param normalized true to learn normalized RBF network.
         */
        public Trainer<T> setNormalized(boolean normalized) {
            this.normalized = normalized;
            return this;
        }

        @Override
        public RBFNetwork<T> train(T[] x, double[] y) {
            @SuppressWarnings("unchecked")
            T[] centers = (T[]) java.lang.reflect.Array.newInstance(x.getClass().getComponentType(), m);
            GaussianRadialBasis gaussian = SmileUtils.learnGaussianRadialBasis(x, centers, distance);
            if (rbf == null) {
                return new RBFNetwork<>(x, y, distance, gaussian, centers, normalized);
            } else {
                return new RBFNetwork<>(x, y, distance, rbf, centers, normalized);
            }
        }
        
        /**
         * Learns a RBF network with given centers.
         * 
         * @param x the training data.
         * @param y the response variable.
         * @param centers the centers of RBF functions.
         * @return a trained RBF network
         */
        public RBFNetwork<T> train(T[] x, double[] y, T[] centers) {
            return new RBFNetwork<>(x, y, distance, rbf, centers, normalized);
        }
    }
    
    /**
     * Constructor. Learn a regular RBF network without normalization.
     * @param x the training data.
     * @param y the response variable.
     * @param distance the distance functor.
     * @param rbf the radial basis function.
     * @param centers the centers of RBF functions.
     */
    public RBFNetwork(T[] x, double[] y, Metric<T> distance, RadialBasisFunction rbf, T[] centers) {
        this(x, y, distance, rbf, centers, false);
    }

    /**
     * Constructor. Learn a regular RBF network without normalization.
     * @param x the training data.
     * @param y the response variable.
     * @param distance the distance functor.
     * @param rbf the radial basis functions.
     * @param centers the centers of RBF functions.
     */
    public RBFNetwork(T[] x, double[] y, Metric<T> distance, RadialBasisFunction[] rbf, T[] centers) {
        this(x, y, distance, rbf, centers, false);
    }

    /**
     * Constructor.
     * @param x the training data.
     * @param y the response variable.
     * @param distance the distance functor.
     * @param rbf the radial basis function.
     * @param centers the centers of RBF functions.
     * @param normalized true for the normalized RBF network.
     */
    public RBFNetwork(T[] x, double[] y, Metric<T> distance, RadialBasisFunction rbf, T[] centers, boolean normalized) {
        this(x, y, distance, rep(rbf, centers.length), centers, normalized);
    }
    
    /**
     * Returns an array of radial basis functions initialized with given values.
     * @param rbf the initial value of array.
     * @param k the size of array.
     * @return an array of radial basis functions initialized with given values
     */
    private static RadialBasisFunction[] rep(RadialBasisFunction rbf, int k) {
        RadialBasisFunction[] arr = new RadialBasisFunction[k];
        Arrays.fill(arr, rbf);
        return arr;
    }
    
    /**
     * Constructor.
     * @param x the training dataset.
     * @param y the response variable.
     * @param distance the distance functor.
     * @param rbf the radial basis functions.
     * @param centers the centers of RBF functions.
     * @param normalized true for the normalized RBF network.
     */
    public RBFNetwork(T[] x, double[] y, Metric<T> distance, RadialBasisFunction[] rbf, T[] centers, boolean normalized) {
        if (x.length != y.length) {
            throw new IllegalArgumentException(String.format("The sizes of X and Y don't match: %d != %d", x.length, y.length));
        }

        if (rbf.length != centers.length) {
            throw new IllegalArgumentException(String.format("The sizes of RBF functions and centers don't match: %d != %d", rbf.length, centers.length));
        }

        this.centers = centers;
        this.distance = distance;
        this.rbf = rbf;
        this.normalized = normalized;
        
        int n = x.length;
        int m = rbf.length;

        DenseMatrix G = Matrix.zeros(n, m);
        double[] b = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < m; j++) {
                double r = rbf[j].f(distance.d(x[i], centers[j]));
                G.set(i, j, r);
                sum += r;
            }

            if (normalized) {
                b[i] = sum * y[i];
            } else {
                b[i] = y[i];
            }
        }

        w = new double[m];
        QR qr = G.qr();
        qr.solve(b, w);
    }

    @Override
    public double predict(T x) {
        double sum = 0.0, sumw = 0.0;
        for (int i = 0; i < rbf.length; i++) {
            double f = rbf[i].f(distance.d(x, centers[i]));//fzhang 2019.8.24 the distance between x and the center
            sumw += w[i] * f;
            sum += f;
        }

        return normalized ? sumw / sum : sumw;
    }

//    public static void main(String[] args){
//        double[][][] x= new double[2][15][3];
//        for(int m = 0; m < 2; m++){
//            for(int i=0; i<15; i++)
//                for(int j=0; j<3; j++)
//                    x[m][i][j] = (i+j+0.0);
//        }
//
//
//        double[] y = new double[]{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};
//
//        Metric<double[]> metric = new EuclideanDistance();
//        Trainer<double[]> trainer = new Trainer<double[]>(metric);
//        RBFNetwork<double[]> network = trainer.train(x[0], y);
//
//        double[] x_test = new double[]{0.0, 1.1, 2.0};
//        System.out.println(network.predict(x_test));


    public static void main(String[] args){
        double[][] x= new double[2][3];
        for(int i=0; i<2; i++)
            for(int j=0; j<3; j++)
                x[i][j] = (i+j+0.0);


        double[] y = new double[]{1.0,2.0};

        Metric<double[]> metric = new EuclideanDistance();
        Trainer<double[]> trainer = new Trainer<double[]>(metric);
        RBFNetwork<double[]> network = trainer.train(x, y);

        double[] x_test = new double[]{0.0, 1.1};
        System.out.println(network.predict(x_test));
    }
}


