/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


/**
 *
 * @author Fei Lu 
 */
public class Bins {
    int corStart;
    /**Exclusive*/
    int corEnd;
    Interval[] bins;
    
    public Bins (int corStart, int corEnd, int rangeLength, int[] coordinate, double[] v) {
        this.initialize(corStart, corEnd, rangeLength);
        this.loadValues(coordinate, v);
    }
    
    public Bins (int corStart, int corEnd, int rangeLength) {
        this.initialize(corStart, corEnd, rangeLength);
    }
    
    private void initialize (int corStart, int corEnd, int rangeLength) {
        this.corStart = corStart;
        this.corEnd = corEnd;
        int maxLength = corEnd-corStart;
        int base = maxLength%rangeLength;
        int binNum;
        if  (base == 0) {
            binNum = maxLength/rangeLength;
        }
        else {
            binNum = maxLength/rangeLength+1;
        }
        bins = new Interval[binNum];
        for (int i = 0; i < binNum; i++) {
            int actualLength = rangeLength;
            if (i == binNum-1 && base != 0) {
                actualLength = base;
            }
            bins[i] = new Interval(i*rangeLength+corStart, i*rangeLength+corStart+actualLength);
        }
        Arrays.sort(bins);
    }
    
    public int getBinNum () {
        return bins.length;
    }
    
    public int getCoordinateStart() {
        return corStart;
    }
    
    public int getCoordinateEnd() {
        return corEnd;
    }
    
    public int getBinStart (int binIndex) {
        return bins[binIndex].start;
    }
    
    public int getBinEnd (int binIndex) {
        return bins[binIndex].end;
    }
    
    public double getBinMedian (int binIndex) {
        return bins[binIndex].mean;
    }
    
    public double getBinAverage (int binIndex) {
        if (!this.isThereValueInBin(binIndex)) return Double.NaN;
        double[] v = this.getBinValues(binIndex);
        double aver = 0;
        for (int i = 0; i < v.length; i++) aver+= v[i]/v.length;
        return aver;
    }
    
    public double getTotalValue (int binIndex) {
        if (!this.isThereValueInBin(binIndex)) return Double.NaN;
        double v = 0;
        for (int i = 0; i < this.getNumValues(binIndex); i++) {
            v+=bins[binIndex].values[i];
        }
        return v;
    }
    
    public int getNumValues (int binIndex) {
        return bins[binIndex].values.length;
    }
    
    public double getBinSD (int binIndex) {
        return bins[binIndex].sd;
    }
    
    public double[] getBinValues (int binIndex) {
        return bins[binIndex].values;
    }
    
    public boolean isThereValueInBin (int binIndex) {
        if (this.getBinValues(binIndex).length == 0) return false;
        return true;
    }
    
    public void loadValues (int[] coordinate, double[] v) {
        ArrayList<Double>[] vs = new ArrayList[bins.length];
        Interval query;
        int index;
        for (int i = 0; i < vs.length; i++) {
            vs[i] = new ArrayList();
        }
        for (int i = 0; i < coordinate.length; i++) {
            query = new Interval(coordinate[i], coordinate[i]);
            index = Arrays.binarySearch(bins, query);
            if (index > 0) {
                
            }
            else if (index < -1 && index >= -bins.length) {
                index = -index-2;
            }
            else if (index == -1) {
                
            }
            else if (index < -bins.length) {
                index = -index-2;
                if (coordinate[i] >= bins[index].end) {
                    index = -1;
                }
            }
            if (!(index < 0)) {
                vs[index].add(v[i]);
            }
        }
        for (int i = 0; i < bins.length; i++) {
            bins[i].loadValues(vs[i]);
        }
    }
    
    public void writeBinSummary (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("BinStart\tBinEnd\tMean\tSD");
            bw.newLine();
            for (int i = 0; i < bins.length; i++) {
                bw.write(String.valueOf(bins[i].start)+"\t"+String.valueOf(bins[i].end)+"\t"+String.valueOf(bins[i].mean)+"\t"+String.valueOf(bins[i].sd));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    class Interval implements Comparable<Interval> {
        int start;
        /**Exclusive*/
        int end;
        double[] values = new double[0];
        double sd = Double.NaN;
        double mean = Double.NaN;
        
        Interval (int start, int end) {
            this.start = start;
            this.end = end;
        }
        
        void loadValues (ArrayList<Double> list) {
            if (list.isEmpty()) {
                values = new double[0];
                sd = Double.NaN;
                mean = Double.NaN;
            }
            else {
                values = new double[list.size()];
                for (int i = 0; i < list.size(); i++) {
                    values[i] = list.get(i);
                }
                DescriptiveStatistics ds = new DescriptiveStatistics(values);
                mean = ds.getMean();
                sd = ds.getStandardDeviation();   
            }
        }
        
        @Override
        public int compareTo(Interval o) {
            return start-o.start;
        }
    }
    
}
