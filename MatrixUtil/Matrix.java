//Version 0.2.0
//Last updated 7/29/2017
//Created by Matthew Ognibene
package MatrixUtil;
//TODO shorts
public final class Matrix {
    private Matrix(){}

    public static double[][] dot(double[][] a, double [][] b){
        if(a[0].length==b.length){
            double[][] newMatrix = new double[a.length][b[0].length];
            for(int rowA = 0; rowA < a.length; rowA++){
                for(int colB = 0; colB < b[0].length; colB++) {
                    double sum = 0;
                    for (int aIterator = 0, bIterator = 0; aIterator < a[0].length; aIterator++, bIterator++) {
                        sum += a[rowA][aIterator] * b[bIterator][colB];
                    }
                    newMatrix[rowA][colB] = sum;
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("The number of columns in the first matrix must be equal to the number of rows in the second matrix");
        }
    }
    public static int[][] dot(int[][] a, int [][] b){
        if(a[0].length==b.length){
            int[][] newMatrix = new int[a.length][b[0].length];
            for(int rowA = 0; rowA < a.length; rowA++){
                for(int colB = 0; colB < b[0].length; colB++) {
                    int sum = 0;
                    for (int aIterator = 0, bIterator = 0; aIterator < a[0].length; aIterator++, bIterator++) {
                        sum += a[rowA][aIterator] * b[bIterator][colB];
                    }
                    newMatrix[rowA][colB] = sum;
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("The number of columns in the first matrix must be equal to the number of rows in the second matrix");
        }
    }
    public static long[][] dot(long[][] a, long [][] b){
        if(a[0].length==b.length){
            long[][] newMatrix = new long[a.length][b[0].length];
            for(int rowA = 0; rowA < a.length; rowA++){
                for(int colB = 0; colB < b[0].length; colB++) {
                    long sum = 0;
                    for (int aIterator = 0, bIterator = 0; aIterator < a[0].length; aIterator++, bIterator++) {
                        sum += a[rowA][aIterator] * b[bIterator][colB];
                    }
                    newMatrix[rowA][colB] = sum;
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("The number of columns in the first matrix must be equal to the number of rows in the second matrix");
        }
    }
    public static float[][] dot(float[][] a, float [][] b){
        if(a[0].length==b.length){
            float[][] newMatrix = new float[a.length][b[0].length];
            for(int rowA = 0; rowA < a.length; rowA++){
                for(int colB = 0; colB < b[0].length; colB++) {
                    float sum = 0;
                    for (int aIterator = 0, bIterator = 0; aIterator < a[0].length; aIterator++, bIterator++) {
                        sum += a[rowA][aIterator] * b[bIterator][colB];
                    }
                    newMatrix[rowA][colB] = sum;
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("The number of columns in the first matrix must be equal to the number of rows in the second matrix");
        }
    }

    //Scalar multiplication
    public static int[][] multiply(int a, int[][] array){
        int[][] newMatrix = new int[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    public static double[][] multiply(double a, int[][] array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    public static double[][] multiply(double a, double[][] array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    public static long[][] multiply(long a, long[][] array){
        long[][] newMatrix = new long[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    public static double[][] multiply(double a, long[][] array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    public static float[][] multiply(float a, float[][] array){
        float[][] newMatrix = new float[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    public static double[][] multiply(double a, float[][] array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }

    public static int[][] add(int[][] a, int [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            int newMatrix[][] = new int[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]+b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }
    public static double[][] add(double[][] a, double [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            double newMatrix[][] = new double[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]+b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }
    public static float[][] add(float[][] a, float [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            float newMatrix[][] = new float[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]+b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }
    public static long[][] add(long[][] a, long [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            long newMatrix[][] = new long[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]+b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }

    public static int[][] subtract(int[][] a, int [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            int newMatrix[][] = new int[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]-b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }
    public static double[][] subtract(double[][] a, double [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            double newMatrix[][] = new double[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]-b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }
    public static long[][] subtract(long[][] a, long [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            long newMatrix[][] = new long[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]-b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }
    public static float[][] subtract(float[][] a, float [][] b){
        if(a.length==b.length&&a[0].length==b[0].length){
            float newMatrix[][] = new float[a.length][a[0].length];
            for(int r = 0; r < a.length; r++) {
                for (int c = 0; c < a[0].length; c++) {
                    newMatrix[r][c] = a[r][c]-b[r][c];
                }
            }
            return newMatrix;
        }else{
            throw new MatrixException("Matrices do not share the same dimensions");
        }
    }

    public static double[][] divide(double[][] a, double[][]b){
        return dot(a,inverse(b));
    }

    public static double[][] inverse(double[][] a){
        //Using minors and cofactors
        if(determinant(a)==0){
            throw new MatrixException("Matrix is singular");
        }
        double[][] matrixOfMinors = new double[a.length][a[0].length];
        for (int r = 0; r < a.length; r++) {
            for (int c = 0; c < a[0].length; c++) {
                double[][]tempMatrix = new double[a.length-1][a.length-1];
                for (int tempR = 0, selectorR = 0; tempR < tempMatrix.length; tempR++, selectorR++) {
                    for (int tempC = 0, selectorC =0; tempC < tempMatrix[0].length; tempC++, selectorC++) {
                        if(selectorR==r){
                            selectorR++;
                        }
                        if(selectorC==c){
                            selectorC++;
                        }
                        tempMatrix[tempR][tempC]=a[selectorR][selectorC];
                    }
                }
                matrixOfMinors[r][c]=determinant(tempMatrix);
            }
        }
        double[][] matrixOfCofactors = new double[matrixOfMinors.length][matrixOfMinors.length];
        for (int r = 0; r < matrixOfCofactors.length; r++) {
            for (int c = 0; c < matrixOfCofactors[0].length; c++) {
                matrixOfCofactors[r][c]=matrixOfMinors[r][c];
                if(r%2==1){
                    matrixOfCofactors[r][c]*=-1;
                }
                if(c%2==1){
                    matrixOfCofactors[r][c]*=-1;
                }
            }
        }
        double[][]adjugate = T(matrixOfCofactors);
        return multiply((1.0/determinant(a)), adjugate);
    }
    public static double[][] inverse(int[][] a){
        return inverse(convertToDoubleMatrix(a));
    }
    public static double[][] inverse(float[][] a){
        return inverse(convertToDoubleMatrix(a));
    }
    public static double[][] inverse(long[][] a){
        return inverse(convertToDoubleMatrix(a));
    }

    //Transpose
    public static int[][] T(int[][]a){
        int[][] newMatrix = new int[a.length][a[0].length];
        for(int r=0, cnew=0; r<a.length;r++,cnew++){
            for(int c=0, rnew=0;c<a[0].length;c++,rnew++){
                newMatrix[rnew][cnew]=a[r][c];
            }
        }
        return newMatrix;
    }
    public static double[][] T(double[][]a){
        double[][] newMatrix = new double[a.length][a[0].length];
        for(int r=0, cnew=0; r<a.length;r++,cnew++){
            for(int c=0, rnew=0;c<a[0].length;c++,rnew++){
                newMatrix[rnew][cnew]=a[r][c];
            }
        }
        return newMatrix;
    }
    public static float[][] T(float[][]a){
        float[][] newMatrix = new float[a.length][a[0].length];
        for(int r=0, cnew=0; r<a.length;r++,cnew++){
            for(int c=0, rnew=0;c<a[0].length;c++,rnew++){
                newMatrix[rnew][cnew]=a[r][c];
            }
        }
        return newMatrix;
    }
    public static long[][] T(long[][]a){
        long[][] newMatrix = new long[a.length][a[0].length];
        for(int r=0, cnew=0; r<a.length;r++,cnew++){
            for(int c=0, rnew=0;c<a[0].length;c++,rnew++){
                newMatrix[rnew][cnew]=a[r][c];
            }
        }
        return newMatrix;
    }

    public static double determinant(double[][] a){
        if(a.length!=a[0].length)
            throw new MatrixException("Matrix is not square");

        if(a.length==2){
            return (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        }else {
            double[] components = new double[a[0].length]; //In the equation (one component per column)
            for (int colIt = 0; colIt < a[0].length; colIt++) {
                double[][] multiplyingMatrix = new double[a.length - 1][a.length - 1];
                for (int r = 0; r < multiplyingMatrix.length; r++) {
                    for (int c = 0, aCol = 0; c < multiplyingMatrix[0].length; c++, aCol++) {
                        if (aCol == colIt) {
                            aCol++;
                        }
                        multiplyingMatrix[r][c] = a[r + 1][aCol];
                    }
                }
                components[colIt] = a[0][colIt] * determinant(multiplyingMatrix);
            }
            for (int i = 0; i < components.length; i++) {
                if (i % 2 == 1) {
                    components[i] *= -1;
                }
            }
            int sum = 0;
            for (double i : components) {
                sum+=i;
            }
            return sum;
        }
    }
    public static double determinant(int[][] a){
        return determinant(convertToDoubleMatrix(a));
    }
    public static double determinant(float[][] a){
        return determinant(convertToDoubleMatrix(a));
    }
    public static double determinant(long[][] a){
        return determinant(convertToDoubleMatrix(a));
    }

    //HERE BE ELEMENTARY ROW OPERATIONS
    public static void switchRows(int row1, int row2, int[][] array){
        int[] temp1 = array[row1];
        array[row1]=array[row2];
        array[row2]=temp1;
    }
    public static void switchRows(int row1, int row2, double[][] array){
        double[] temp1 = array[row1];
        array[row1]=array[row2];
        array[row2]=temp1;
    }
    public static void switchRows(int row1, int row2, float[][] array){
        float[] temp1 = array[row1];
        array[row1]=array[row2];
        array[row2]=temp1;
    }
    public static void switchRows(int row1, int row2, long[][] array){
        long[] temp1 = array[row1];
        array[row1]=array[row2];
        array[row2]=temp1;
    }

    public static void multiplyRow(int row, double arg, double[][] array){
        for(int c = 0; c < array[row].length; c++){
            array[row][c] = array[row][c]*arg;
        }
    }
    public static void multiplyRow(int row, int arg, int[][] array){
        for(int c = 0; c < array[row].length; c++){
            array[row][c] = array[row][c]*arg;
        }
    }
    public static void multiplyRow(int row, float arg, float[][] array){
        for(int c = 0; c < array[row].length; c++){
            array[row][c] = array[row][c]*arg;
        }
    }
    public static void multiplyRow(int row, long arg, long[][] array){
        for(int c = 0; c < array[row].length; c++){
            array[row][c] = array[row][c]*arg;
        }
    }

    public static void addRows(int row1, int row2, int[][]array){
        int[]sum = new int[array.length];
        for(int c = 0 ; c<sum.length; c++){
            sum[c]=array[row1][c]+array[row2][c];
        }
        array[row2]=sum;
    }
    public static void addRows(int row1, int row2, double[][]array){
        double[]sum = new double[array.length];
        for(int c = 0 ; c<sum.length; c++){
            sum[c]=array[row1][c]+array[row2][c];
        }
        array[row2]=sum;
    }
    public static void addRows(int row1, int row2, float[][]array){
        float[]sum = new float[array.length];
        for(int c = 0 ; c<sum.length; c++){
            sum[c]=array[row1][c]+array[row2][c];
        }
        array[row2]=sum;
    }
    public static void addRows(int row1, int row2, long[][]array){
        long[]sum = new long[array.length];
        for(int c = 0 ; c<sum.length; c++){
            sum[c]=array[row1][c]+array[row2][c];
        }
        array[row2]=sum;
    }

    public static double[][] convertToDoubleMatrix(int[][]array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = (double) array[r][c];
            }
        }
        return newMatrix;
    }
    public static double[][] convertToDoubleMatrix(float[][]array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = (double) array[r][c];
            }
        }
        return newMatrix;
    }
    public static double[][] convertToDoubleMatrix(long[][]array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = (double) array[r][c];
            }
        }
        return newMatrix;
    }

    public static void printMatrix(int[][] array){
        for(int r = 0; r < array.length; r++){
            for (int c = 0; c < array[0].length; c++){
                System.out.print(array[r][c] + "\t");
            }
            System.out.println();
        }
    }
    public static void printMatrix(double[][] array){
        for(int r = 0; r < array.length; r++){
            for (int c = 0; c < array[0].length; c++){
                System.out.print(array[r][c] + "\t");
            }
            System.out.println();
        }
    }
    public static void printMatrix(float[][] array){
        for(int r = 0; r < array.length; r++){
            for (int c = 0; c < array[0].length; c++){
                System.out.print(array[r][c] + "\t");
            }
            System.out.println();
        }
    }
    public static void printMatrix(long[][] array){
        for(int r = 0; r < array.length; r++){
            for (int c = 0; c < array[0].length; c++){
                System.out.print(array[r][c] + "\t");
            }
            System.out.println();
        }
    }

    private static class MatrixException extends RuntimeException{
        public MatrixException() {}
        public MatrixException (String message) {super (message);}
    }
}
