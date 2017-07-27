//Last updated 7/27/2017
//Created by Matthew Ognibene
package MatrixUtil;

public final class Matrix {
    private Matrix(){}

    public static double[][] multiply(double[][] a, double [][] b){
        if(a[0].length==b.length){
            double[][] newMatrix = new double[a.length][b[0].length];//TODO
            for(int rowA = 0; rowA < a.length; rowA++){
                for(int colB = 0; colB < b[0].length; colB++) {
                    double sum = 0;//TODO
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
    /*TODO
    public static double[][] multiply(int[][] a, int [][] b){
        return multiply(a,b);
    }
    */

    //Scalar multiplication
    public static double[][] multiply(double a, double[][] array){
        double[][] newMatrix = new double[array.length][array[0].length];
        for(int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[0].length; c++) {
                newMatrix[r][c] = array[r][c]*a;
            }
        }
        return newMatrix;
    }
    /*TODO
    public static int[][] multiply(int a, int[][] array){
        return multiply(a,array);
    }
    public static double[][] multiply(double a, double[][] array){
        return multiply(a,array);
    }
    */
    
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

    public static double[][] divide(double[][] a, double[][]b){
        return multiply(a,inverse(b));
    }

    public static double[][] inverse(double[][] a){//Using minors and cofactors
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
    public static double[][] T(double[][]a){
        double[][] newMatrix = new double[a.length][a[0].length];
        for(int r=0, cnew=0; r<a.length;r++,cnew++){
            for(int c=0, rnew=0;c<a[0].length;c++,rnew++){
                newMatrix[rnew][cnew]=a[r][c];
            }
        }
        return newMatrix;
    }
    /*TODO
    public static double[][] T(int[][]a){
        return T(a);
    }
    */

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
                    for (int c = 0, aCol = 0; c < multiplyingMatrix[0].length; c++, aCol++) {//todo better variable
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

    /*TODO
    public static double determinant(int[][] a){
        return determinant(a);
    }
    */

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

    private static class MatrixException extends RuntimeException{
        public MatrixException() {}
        public MatrixException (String message) {super (message);}
    }
}
