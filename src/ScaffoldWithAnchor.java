/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/


import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;


/**
 *
 * @author Fei Lu
 */
public class ScaffoldWithAnchor {
    ScaffoldAnchorInfo[] sai = null;
    
    /**
     * Instructor from SimpleAlignment and TagsOnGeneticMap, query in SimpleAlignment must be the index of TOGM
     * @param sa
     * @param togm
     */
    public ScaffoldWithAnchor (ShortreadAlignment sa, TagsOnGeneticMap togm, String scaffoldFileS, boolean ifOnlyPerfectMatch) {
        Fasta f = new Fasta (scaffoldFileS);
        f.sortRecordByName();
        sai = new ScaffoldAnchorInfo[f.getSeqNumber()];
        String[] scaffoldNames = f.getNames();
        ArrayList<Integer>[] indicesList = new ArrayList[f.getSeqNumber()];
        for (int i = 0; i < indicesList.length; i++) indicesList[i] = new ArrayList();
        
        ArrayList<String> onlyPerfectMatchList = new ArrayList();
        String[] onlyPerfectMatch;
        if (ifOnlyPerfectMatch) {
            sa.sortByQuery();
            String[] queries = sa.getQuerys();
            for (int i = 0; i < queries.length; i++) {
                if (!sa.isOnlyPerfectMatch(queries[i])) continue;
                onlyPerfectMatchList.add(queries[i]);
            }
            onlyPerfectMatch = onlyPerfectMatchList.toArray(new String[onlyPerfectMatchList.size()]);
            Arrays.sort(onlyPerfectMatch);
            sa.sortByHitAndPos();
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                if (!sa.isPerfectMatch(i)) continue;
                if (Arrays.binarySearch(onlyPerfectMatch, sa.getQuery(i)) < 0) continue;
                int index = Arrays.binarySearch(scaffoldNames, sa.getHit(i));
                indicesList[index].add(i);
            }
        }
        else {
            sa.sortByHitAndPos();
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                if (!sa.isPerfectMatch(i)) continue;
                int index = Arrays.binarySearch(scaffoldNames, sa.getHit(i));
                indicesList[index].add(i);
            }
        }
        for (int i = 0; i < sai.length; i++) {
            Integer[] indices = indicesList[i].toArray(new Integer[indicesList[i].size()]);
            int[] anchorIndex = new int[indices.length];
            int[] anchorGChr = new int[indices.length];
            int[] anchorGPos = new int[indices.length];
            int[] anchorSPPos = new int[indices.length];
            byte[] ifPAV = new byte[indices.length];
            byte[] anchorReadLength = new byte[indices.length];
            for (int j = 0; j < indices.length; j++) {
                anchorIndex[j] = Integer.valueOf(sa.getQuery(indices[j]));
                anchorGChr[j] = togm.getGChr(anchorIndex[j]);
                anchorGPos[j] = togm.getGPos(anchorIndex[j]);
                anchorSPPos[j] = sa.getStartPos(indices[j]);
                ifPAV[j] = togm.getIfPAV(anchorIndex[j]);
                anchorReadLength[j] = (byte)togm.getTagLength(anchorIndex[j]);
            }
            sai[i] = new ScaffoldAnchorInfo(f.getName(i), f.getSeqLength(i), indices.length, anchorIndex, anchorGChr, anchorGPos, anchorSPPos, ifPAV, anchorReadLength);
        }
        Arrays.sort(sai);
    }
    
    public ScaffoldWithAnchor (String scaffoldWithAnchorFileS) {
        this.readScaffoldWithAnchor(scaffoldWithAnchorFileS);
    }
    
    /**
     * Return number of scaffolds
     * @return
     */
    public int getScaffoldNumber () {
        return sai.length;
    }
    
    /**
     * Write scaffold quality
     * @param outfileS
     */
    public void writeScaffoldQuality (String outfileS, String pdfDirS, int minAnchorNum, int[] chromLength, int[] centStart, int[] centEnd) {
        int baseBinSize = 5000000;
        double minProportionAnchorInBaseBin = 0.10;
        double minProprotionAnchorInRegion = 0.2;
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Index\tScaffold\tLength\tAnchorNumber\tMainChromosome\tNumAnchorOnMainChr\tRatioAnchorOnMainChr\tNumRegions\tRegionChr-RegionPosition-RatioAnchorOnRegion");
            bw.newLine();
            Bins[] chromBins = new Bins[chromLength.length];
            ArrayList<Integer>[] chromAnchorList = new ArrayList[chromLength.length];
            for (int i = 0; i < chromBins.length; i++) {
                chromBins[i] = new Bins(1, chromLength[i], baseBinSize);
            }
            int cntIndex = 0;
            for (int i = 0; i < this.getScaffoldNumber(); i++) {
                ScaffoldAnchorInfo ai = sai[i] ;
                if (ai.anchorNum < minAnchorNum)  continue;
                
                String pdfFileS = FStringUtils.getNDigitNumber(7, cntIndex)+"_"+ai.scaffoldName+".pdf";
                pdfFileS = new File(pdfDirS,pdfFileS).getAbsolutePath();
                this.drawPDF(ai, chromLength, centStart, pdfFileS);
                
                int[] chrIndexAnchorCnt = new int[chromLength.length];
                for (int j = 0; j < ai.anchorNum; j++) {
                    chrIndexAnchorCnt[ai.anchorGChr[j]-1]++;
                }
                int mostChrIndex = 0;
                int most = chrIndexAnchorCnt[0];
                for (int j = 1; j < chrIndexAnchorCnt.length; j++) {
                    if (chrIndexAnchorCnt[j] > most) {
                        most = chrIndexAnchorCnt[j];
                        mostChrIndex = j;
                    }
                }
                int mainChr = mostChrIndex+1;
                double ratioAnchorOnMainChr= (double)chrIndexAnchorCnt[mostChrIndex]/ai.anchorNum;
                for (int j = 0; j < chromAnchorList.length; j++) {
                    chromAnchorList[j] = new ArrayList();
                }
                for (int j = 0; j < ai.anchorNum; j++) {
                    chromAnchorList[ai.anchorGChr[j]-1].add(ai.anchorGPos[j]);
                }
                for (int j = 0; j < chromAnchorList.length; j++) {
                    int[] coor = new int[chromAnchorList[j].size()];
                    double[] value = new double[coor.length];
                    for (int k = 0; k < coor.length; k++) {
                        coor[k] = chromAnchorList[j].get(k);
                        value[k] = 1;
                    }
                    chromBins[j].loadValues(coor, value);
                }
                int anchorCut = (int)(ai.anchorNum*minProportionAnchorInBaseBin);
                ArrayList<Integer> passChrIndexList = new ArrayList();
                ArrayList<Integer> passBinIndexList = new ArrayList();
                ArrayList<Integer> passBinAnchorNum = new ArrayList();
                for (int j = 0; j < chromBins.length; j++) {
                    for (int k = 0; k < chromBins[j].getBinNum(); k++) {
                        int nu = chromBins[j].getNumValues(k);
                        if (nu == 0) continue;
                        if (nu < anchorCut) continue;
                        passChrIndexList.add(j);
                        passBinIndexList.add(k);
                        passBinAnchorNum.add(nu);
                    }
                }
                ArrayList<ArrayList<Integer>> regionIndicesList = new ArrayList();
                
                int cnt = 0;
                
                if (!passBinIndexList.isEmpty()) {
                    regionIndicesList.add(new ArrayList());
                    regionIndicesList.get(cnt).add(0);
                }
                for (int j = 0; j < passBinIndexList.size()-1; j++) {
                    if (((passBinIndexList.get(j)-passBinIndexList.get(j+1))) == -1 && passChrIndexList.get(j) == passChrIndexList.get(j+1)) {
                        regionIndicesList.get(cnt).add(j+1);
                    }
                    else {
                        cnt++;
                        regionIndicesList.add(new ArrayList());
                        regionIndicesList.get(cnt).add(j+1);
                    }
                }
                int numRegions = regionIndicesList.size();
                int[] anchorNumInRegion = new int[numRegions];
                int[] regionChr = new int[numRegions];
                int[] regionPos = new int[numRegions];
                double[] ratioAnchorOnRegion = new double[numRegions];
                for (int j = 0; j < numRegions; j++) {
                    Integer[] regionIndices = regionIndicesList.get(j).toArray(new Integer[regionIndicesList.get(j).size()]);
                    int chrIndex = 0;
                    try {
                        chrIndex = passChrIndexList.get(regionIndicesList.get(j).get(0));
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                    regionChr[j] = chrIndex+1;
                    int[] binIndex = new int[regionIndices.length];
                    for (int k = 0; k < binIndex.length; k++) {
                        binIndex[k] = passBinIndexList.get(regionIndicesList.get(j).get(k));
                    }
                    regionPos[j] = (binIndex[0]+1+binIndex[binIndex.length-1]) * baseBinSize/2;
                    for (int k = 0; k < binIndex.length; k++) {
                        anchorNumInRegion[j]+= chromBins[chrIndex].getNumValues(binIndex[k]);
                    }
                    if (binIndex[0] > 0) anchorNumInRegion[j] += chromBins[chrIndex].getNumValues(binIndex[0]-1);
                    if (binIndex[binIndex.length-1] < chromBins[chrIndex].getBinNum()-1) anchorNumInRegion[j] += chromBins[chrIndex].getNumValues(binIndex[binIndex.length-1]+1);
                    ratioAnchorOnRegion[j] = (double)anchorNumInRegion[j]/ai.anchorNum;
                }
                bw.write(String.valueOf(cntIndex)+"\t"+ai.scaffoldName+"\t"+String.valueOf(ai.scaffoldLength)+"\t"+String.valueOf(ai.anchorNum)+"\t"+String.valueOf(mainChr)+"\t"+String.valueOf(chrIndexAnchorCnt[mostChrIndex])+"\t"+get3DigitDouble(ratioAnchorOnMainChr)+"\t");
                bw.write(String.valueOf(numRegions)+"\t");
                for (int j = 0; j < numRegions-1; j++) {
                    bw.write(String.valueOf(regionChr[j])+"-"+String.valueOf(regionPos[j])+"-"+get3DigitDouble(ratioAnchorOnRegion[j])+";");
                }
                if (numRegions>0) {
                    bw.write(String.valueOf(regionChr[numRegions-1])+"-"+String.valueOf(regionPos[numRegions-1])+"-"+String.valueOf(anchorNumInRegion[numRegions-1])+"-"+get3DigitDouble(ratioAnchorOnRegion[numRegions-1]));
                }
                bw.newLine();
                cntIndex++;
            }
            bw.flush();
            bw.close();
            System.out.println("ScaffoldQuality file is written to "  + outfileS);
            System.out.println("Evaluation figures are in "  + pdfDirS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Write scaffold quality
     * @param outfileS
     */
    public void writeScaffoldQualitySimple (String outfileS, String pdfDirS, int minAnchorNum, int[] chromLength, int[] centStart, int[] centEnd) {

        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Index\tScaffold\tLength\tAnchorNumber\tMainChromosome\tNumAnchorOnMainChr\tRatioAnchorOnMainChr");
            bw.newLine();

            ArrayList<Integer>[] chromAnchorList = new ArrayList[chromLength.length];
            
            int cntIndex = 0;
            for (int i = 0; i < this.getScaffoldNumber(); i++) {
                ScaffoldAnchorInfo ai = sai[i] ;
                if (ai.anchorNum < minAnchorNum)  continue;
                
                String pdfFileS = FStringUtils.getNDigitNumber(7, cntIndex)+"_"+ai.scaffoldName+".pdf";
                pdfFileS = new File(pdfDirS,pdfFileS).getAbsolutePath();
                this.drawPDF(ai, chromLength, centStart, pdfFileS);
                
                int[] chrIndexAnchorCnt = new int[chromLength.length];
                for (int j = 0; j < ai.anchorNum; j++) {
                    chrIndexAnchorCnt[ai.anchorGChr[j]-1]++;
                }
                int mostChrIndex = 0;
                int most = chrIndexAnchorCnt[0];
                for (int j = 1; j < chrIndexAnchorCnt.length; j++) {
                    if (chrIndexAnchorCnt[j] > most) {
                        most = chrIndexAnchorCnt[j];
                        mostChrIndex = j;
                    }
                }
                int mainChr = mostChrIndex+1;
                double ratioAnchorOnMainChr= (double)chrIndexAnchorCnt[mostChrIndex]/ai.anchorNum;
                for (int j = 0; j < chromAnchorList.length; j++) {
                    chromAnchorList[j] = new ArrayList();
                }
                for (int j = 0; j < ai.anchorNum; j++) {
                    chromAnchorList[ai.anchorGChr[j]-1].add(ai.anchorGPos[j]);
                } 
                bw.write(String.valueOf(cntIndex)+"\t"+ai.scaffoldName+"\t"+String.valueOf(ai.scaffoldLength)+"\t"+String.valueOf(ai.anchorNum)+"\t"+String.valueOf(mainChr)+"\t"+String.valueOf(chrIndexAnchorCnt[mostChrIndex])+"\t"+get3DigitDouble(ratioAnchorOnMainChr));
                bw.newLine();
                cntIndex++;
            }
            bw.flush();
            bw.close();
            System.out.println("ScaffoldQuality file is written to "  + outfileS);
            System.out.println("Evaluation figures are in "  + pdfDirS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private String get3DigitDouble (double value) {
        DecimalFormat df = new DecimalFormat("#.000");
        return "0"+df.format(value);
    }
    
    private void drawPDF (ScaffoldAnchorInfo ai, int[] chromLength, int[] centromerePos, String pdfFileS) {
        int chrNum = chromLength.length;
        int[] chrID = new int[chrNum];
        for (int i = 0; i < chrID.length; i++) chrID[i] = i+1;
        String scaffoldName = ai.scaffoldName;
        int scaffoldLength = ai.scaffoldLength;
        int anchorNum = ai.anchorNum;
        int[] anchorScaffoldPos = Arrays.copyOf(ai.anchorSPPos, ai.anchorNum);
        int[] anchorGChr = Arrays.copyOf(ai.anchorGChr, ai.anchorNum);
        int[] anchorGPos = Arrays.copyOf(ai.anchorGPos, ai.anchorNum);
        int[] chrLength = Arrays.copyOf(chromLength, chrNum);
        int[] centPos = Arrays.copyOf(centromerePos, chrNum);
        
        int width = 1200;
        int height = 400;
        int maxRecLength = 200;
        int xStart = 50;
        int xInterval = 20;
        int recHeight = 6;
        int yStart = 50;
        int yInterval = 150;
        int contigRecLength = maxRecLength*5+xInterval*4;
        double ratio = (double)maxRecLength/chrLength[0];
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = (int)(ratio*chrLength[i]);
            centPos[i] = (int)(ratio*centPos[i]);
        }
        double ratioScaffold = (double)contigRecLength/scaffoldLength;
        for (int i = 0; i < anchorNum; i++) {
            anchorScaffoldPos[i] = (int)(anchorScaffoldPos[i]*ratioScaffold)+xStart;
            int n = (anchorGChr[i]-1)%5;
            anchorGPos[i] = xStart+n*(maxRecLength+xInterval)+(int)(anchorGPos[i]*ratio);
        }
        try {
            Document pd = new Document(new Rectangle(width,height));
            PdfWriter pw;
            pw = PdfWriter.getInstance(pd, new FileOutputStream (pdfFileS));
            pd.open();
            PdfContentByte canvas = pw.getDirectContent();
            DefaultFontMapper mapper = new DefaultFontMapper();
            PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
            int x = xStart;
            int y = yStart;
            g2d.drawString(scaffoldName+"\t\tLength="+String.valueOf(scaffoldLength), x, y-25);
            for (int j = 0; j < 5; j++) {
                g2d.setColor(Color.black);
                g2d.drawString(String.valueOf(j+1), x, y-5);
                g2d.setColor(Color.ORANGE);
                g2d.fillRect(x, y, chrLength[j], recHeight);
                g2d.setColor(Color.black);
                g2d.fillRect(x+centPos[j], y, 3, recHeight);
                x+=chrLength[0];
                x+=xInterval;
            }
            x = xStart;
            y+=recHeight; int y1 = y;
            y+=yInterval; int y2 = y;
            g2d.setColor(Color.cyan);
            g2d.fillRect(x, y, contigRecLength, recHeight);
            y+=recHeight; int y3 = y;
            y+=yInterval; int y4 = y;
            for (int j = 5; j < chrNum; j++) {
                g2d.setColor(Color.black);
                g2d.drawString(String.valueOf(j+1), x, y+recHeight+15);
                g2d.setColor(Color.ORANGE);
                g2d.fillRect(x, y, chrLength[j], recHeight);
                g2d.setColor(Color.black);
                g2d.fillRect(x+centPos[j], y, 3, recHeight);
                x+=chrLength[0];
                x+=xInterval;
            }
            Color tBlue = new Color(0, 0, 255, 50);
            g2d.setColor(tBlue);
            for (int i = 0; i < anchorNum; i++) {
                if (anchorGChr[i] <= 5){
                    g2d.drawLine(anchorGPos[i], y1, anchorScaffoldPos[i], y2);
                }
                else {
                    g2d.drawLine(anchorScaffoldPos[i], y3,anchorGPos[i] , y4);
                }
            }
            g2d.dispose();
            pd.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    /**
     * Write scaffoldWithAnchor
     * @param outfileS
     */
    public void writeScaffoldWithAnchor (String outfileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("#>Scaffold\tScaffoldLength\tAnchorNumber");
            bw.newLine();
            bw.write("#AnchorIndex\tScaffoldPosition\tGeneticChr\tGeneticPosition\tIfPAV\tAnchorReadLength");
            bw.newLine();
            for (int i = 0; i < sai.length; i++) {
                bw.write(">"+sai[i].scaffoldName+"\t"+String.valueOf(sai[i].scaffoldLength)+"\t"+String.valueOf(sai[i].anchorNum));
                bw.newLine();
                for (int j = 0; j < sai[i].anchorNum; j++) {
                    bw.write(String.valueOf(sai[i].anchorIndex[j])+"\t"+String.valueOf(sai[i].anchorSPPos[j])+"\t"+sai[i].anchorGChr[j]+"\t"+String.valueOf(sai[i].anchorGPos[j]));
                    bw.write("\t"+String.valueOf(sai[i].ifPAV[j])+"\t"+String.valueOf(sai[i].anchorReadLength[j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("ScaffoldWithAnchor file written to " + outfileS);
    }
    
    public void readScaffoldWithAnchor (String infileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            int cnt = 0;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) cnt++;
            }
            sai = new ScaffoldAnchorInfo[cnt];
            br = IoUtils.getTextReader(infileS);
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    temp = temp.replaceFirst(">", "");
                    String[] tem = temp.split("\t");
                    String scaffoldName = tem[0];
                    int scaffoldLength = Integer.valueOf(tem[1]);
                    int anchorNum = Integer.valueOf(tem[2]);
                    int[] anchorIndex = new int[anchorNum];
                    int[] anchorGChr = new int[anchorNum];
                    int[] anchorGPos = new int[anchorNum];
                    int[] anchorSPPos = new int[anchorNum];
                    byte[] ifPAV = new byte[anchorNum];
                    byte[] anchorReadLength = new byte[anchorNum];
                    for (int i = 0; i < anchorNum; i++) {
                        tem = br.readLine().split("\t");
                        anchorIndex[i] = Integer.valueOf(tem[0]);
                        anchorSPPos[i] = Integer.valueOf(tem[1]);
                        anchorGChr[i] = Integer.valueOf(tem[2]);
                        anchorGPos[i] = Integer.valueOf(tem[3]);
                        ifPAV[i] = Byte.valueOf(tem[4]);
                        anchorReadLength[i] = Byte.valueOf(tem[5]);
                    }
                    sai[cnt] = new ScaffoldAnchorInfo(scaffoldName, scaffoldLength, anchorNum, anchorIndex, anchorGChr, anchorGPos, anchorSPPos, ifPAV, anchorReadLength);
                    cnt++;
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("ScaffoldWithAnchor file read from " + infileS);
    }
    
    class ScaffoldAnchorInfo implements Comparable<ScaffoldAnchorInfo> {
        String scaffoldName = null;
        int scaffoldLength = Integer.MIN_VALUE;
        int anchorNum = Integer.MIN_VALUE;
        int[] anchorIndex = null;
        /**Physical position of anchors on scaffold*/
        int[] anchorSPPos = null;
        int[] anchorGChr = null;
        int[] anchorGPos = null;
        byte[] ifPAV = null;
        byte[] anchorReadLength = null;
        
        ScaffoldAnchorInfo (String scaffoldName, int scaffoldLength, int anchorNum, int[] anchorIndex, int[] anchorGChr, int[] anchorGPos, int[] anchorSPPos, byte[] ifPAV, byte[] anchorReadLength) {
            this.scaffoldName = scaffoldName;
            this.scaffoldLength = scaffoldLength;
            this.anchorNum = anchorNum;
            this.anchorIndex = anchorIndex;
            this.anchorGChr = anchorGChr;
            this.anchorGPos = anchorGPos;
            this.anchorSPPos = anchorSPPos;
            this.ifPAV = ifPAV;
            this.anchorReadLength = anchorReadLength;
        }
        
        @Override
        public int compareTo(ScaffoldAnchorInfo o) {
            if (this.scaffoldLength > o.scaffoldLength) return -1;
            else if (this.scaffoldLength < o.scaffoldLength) return 1;
            else return 0;
        }
    }
    
}
