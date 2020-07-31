package org.ohdsi.data;

/**
 * @author Marc A. Suchard
 */
public class SortedCoxData {

    public int[] y;
    public double[] x;
    public int[] strata;
    public double[] weight;

    public SortedCoxData(int[] y, double[] x, int[] strata, double[] weight) {
        this.y = y;
        this.x = x;
        this.strata = strata;
        this.weight = weight;
    }
}
