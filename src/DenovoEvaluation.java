/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/




import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;


/**
 * This pipeline is used to evaluate de novo assembly using pan-genome anchors. The input is only a parameter file. Run on linux
 * @author Fei Lu
 */
public class DenovoEvaluation {
    String workingDirS = null;
    String assemblySourceFileS = null;
    String togmFileS = null;
    int minScaffoldLength = 20000;
    int minAnchorNum = 20;
    
    String[] subDirS = {"genome", "alignment", "pdf", "result"};
    BufferedWriter bwLog = null;
    String genomeDirS = null;
    String alignmentDirS = null;
    String pdfDirS = null;
    String resultDirS = null;
    String rPath = null;
    
    int[] chromLength =null;
    int[] centStart = null;
    int[] centEnd = null;
    
    public DenovoEvaluation (String parameterFileS) {
        this.getParameters(parameterFileS);
        this.pipe();
    }
    
    public DenovoEvaluation () {
        String parameterFileS = "M:\\pipelineTest\\denovoEvaluation\\smallset\\parameter.txt";
        this.getParameters(parameterFileS);
        this.test();
    }
    
    private void pipe () {
        this.mkDirS();
        this.selectScaffold();
        this.togmToFastq();
        this.bowtie2Build();
        this.alignBowtie2();
        this.convertAlignment();
        this.outputResult();
        System.out.println("Denovo assembly evaluation finished");
    }
    
    private void test () {
        String scaffoldWithAnchorFileS = "M:\\pipelineTest\\denovoEvaluation\\smallset\\cml247Evaluation\\result\\scaffoldWithAnchor.txt";
        ScaffoldWithAnchor swa = new ScaffoldWithAnchor(scaffoldWithAnchorFileS);
        String outfileName = "scaffoldQuality_"+String.valueOf(minScaffoldLength)+"bp_"+String.valueOf(this.minAnchorNum)+"anchors.txt";
        String scaffoldQualityFileS = new File("M:\\pipelineTest\\denovoEvaluation\\smallset\\cml247Evaluation\\result\\", outfileName).getAbsolutePath();
        String pdfDirS = "M:\\pipelineTest\\denovoEvaluation\\smallset\\cml247Evaluation\\pdf\\";
        swa.writeScaffoldQualitySimple(scaffoldQualityFileS, pdfDirS, minAnchorNum, this.chromLength, this.centStart, this.centEnd);
    }
    
    private void outputResult () {
        String selectedScaffoldsFileS = new File(genomeDirS, "selectedScaffolds.fas").getAbsolutePath();
        String alignmentFileS = new File(alignmentDirS, "alignment.sa.bin").getAbsolutePath();
        String scaffoldWithAnchorFileS = new File(resultDirS, "scaffoldWithAnchor.txt").getAbsolutePath();
        ShortreadAlignment sa = new ShortreadAlignment(alignmentFileS, IOFileFormat.Binary);
        TagsOnGeneticMap togm = new TagsOnGeneticMap(togmFileS, TagsOnGeneticMap.FilePacking.Text);
        ScaffoldWithAnchor swa = new ScaffoldWithAnchor(sa, togm, selectedScaffoldsFileS, true);
        swa.writeScaffoldWithAnchor(scaffoldWithAnchorFileS);
        this.outputLog("ScaffoldWithAnchor file is written to "  + scaffoldWithAnchorFileS + "\n");
        String outfileName = "scaffoldQuality_"+String.valueOf(minScaffoldLength)+"bp_"+String.valueOf(this.minAnchorNum)+"anchors.txt";
        String scaffoldQualityFileS = new File(resultDirS, outfileName).getAbsolutePath();
        swa.writeScaffoldQualitySimple(scaffoldQualityFileS, pdfDirS, minAnchorNum, this.chromLength, this.centStart, this.centEnd);
        
        String scaffoldQualityPdf = scaffoldQualityFileS.replaceFirst(".txt", ".pdf");
        Table t = new Table (scaffoldQualityFileS);
        double[] v = new double[t.getRowNumber()];
        for (int i = 0; i < v.length; i++) v[i] = Double.valueOf(t.content[i][6]);
        Histogram d = new Histogram(v);
        d.setRPath(rPath);
        d.setTitle("Assembly quality");
        d.setXLab("RatioAnchorOnMainChr");
        d.setYLab("Frequency");
        d.setSlideMode();
        d.saveGraph(scaffoldQualityPdf);
        this.outputLog("ScaffoldQuality file is written to "  + scaffoldQualityFileS + "\n");
        this.outputLog("Evaluation figures are in "  + pdfDirS + "\n\n");
    }
    
    private void convertAlignment () {
        String samFileS = new File(alignmentDirS, "alignmentK2.sam").getAbsolutePath();
        String alignmentFileS = new File(alignmentDirS, "alignment.sa.bin").getAbsolutePath();
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        sa.writeSimpleAlignment(alignmentFileS, IOFileFormat.Binary);
        System.out.println("SAM file is converted to " + alignmentFileS);
        this.outputLog("SAM file is converted to " + alignmentFileS + "\n\n");
    }
    
    private void alignBowtie2 () {
        try {
            String libFileS = new File(genomeDirS, "selectedScaffolds.fas").getAbsolutePath();
            String queryFileS = new File(alignmentDirS, "anchor.fq").getAbsolutePath();
            String samFileS = new File(alignmentDirS, "alignmentK2.sam").getAbsolutePath();
            int numOfCores = Runtime.getRuntime().availableProcessors();
            String cmd = "bowtie2 -x " + libFileS + " -q " + queryFileS + " -k 2 --very-sensitive-local -S " + samFileS + " -p " + String.valueOf(numOfCores);
            Process p = Runtime.getRuntime().exec(cmd);
            System.out.println(cmd);
            System.out.println("Alignment using Bowtie2. Please wait...");
            System.out.println("Patience is bitter but its fruit is sweet...");
            p.waitFor();
            System.out.println("Alignment finished");
            this.outputLog(cmd+"\n");
            this.outputLog("SAM file output to " + samFileS + "\n\n");
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void bowtie2Build () {
        try {
            String inputFileS = new File(genomeDirS, "selectedScaffolds.fas").getAbsolutePath();
            String cmd = "bowtie2-build " + inputFileS + " " + inputFileS;
            Process p = Runtime.getRuntime().exec(cmd);
            System.out.println("Indexing scaffolds using Bowtie2...");
            p.waitFor();
            System.out.println("Indexing finished");
            this.outputLog("Indexing finished at " + inputFileS + "\n\n");
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void togmToFastq () {
        String anchorFastqFileS = new File(alignmentDirS, "anchor.fq").getAbsolutePath();
        this.outputLog("TOGM file is at " + anchorFastqFileS + "\n");
        TagsOnGeneticMap togm = new TagsOnGeneticMap (togmFileS, TagsOnGeneticMap.FilePacking.Text);
        togm.writeFastQ(anchorFastqFileS);
        this.outputLog("TOGM file is converted to FASTQ at " + anchorFastqFileS + "\n\n");
    }
    
    private void selectScaffold () {
        this.outputLog("Reading genome assembly file from " + assemblySourceFileS+ "\n\n");
        Fasta f = new Fasta (assemblySourceFileS);
        f.sortRecordByLengthDescending();
        this.stackLog("Genome statistics:\n");
        this.stackLog("Total length: " + String.valueOf(f.getTotalSeqLength()) + "\n");
        this.stackLog("Total scaffold number: " + String.valueOf(f.getSeqNumber())+"\n");
        this.stackLog("L50 length: " + String.valueOf(f.getL50()) + " bp\n");
        this.stackLog("N50 number: " + String.valueOf(f.getN50()) + "\n");
        this.flushLog();
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        for (int i = 0; i < ifOut.length; i++) {
            if (f.getSeqLength(i) > this.minScaffoldLength) ifOut[i] = true;
        }
        String outfileS = new File(genomeDirS, "selectedScaffolds.fas").getAbsolutePath();
        f.writeFasta(outfileS, ifOut);
        this.outputLog("Selected scaffolds are output to " + outfileS +"\n\n");
    }
    
    private void mkDirS () {
        File workingDir = new File(workingDirS);
        workingDir.mkdir();
        for (int i = 0; i < this.subDirS.length; i++) {
            new File(workingDirS, subDirS[i]).mkdir();
        }
        try {
            bwLog = IoUtils.getTextWriter(new File(workingDir, "result_log.txt").getAbsolutePath());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        genomeDirS = new File(workingDirS, subDirS[0]).getAbsolutePath();
        alignmentDirS = new File(workingDirS, subDirS[1]).getAbsolutePath();
        pdfDirS = new File(workingDirS, subDirS[2]).getAbsolutePath();
        resultDirS = new File(workingDirS, subDirS[3]).getAbsolutePath();
    }
    
    private void outputLog (String s) {
        stackLog(s);
        flushLog();
    }
    
    private void stackLog (String s) {
        try {
            bwLog.write(s);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void flushLog () {
        try {
            bwLog.flush();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void getParameters (String parameterFileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(parameterFileS);
            String temp;
            ArrayList<String> pList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                pList.add(temp);
            }
            workingDirS = pList.get(0);
            assemblySourceFileS = pList.get(1);
            togmFileS = pList.get(2);
            this.minScaffoldLength = Integer.valueOf(pList.get(3));
            this.minAnchorNum = Integer.valueOf(pList.get(4));
            this.rPath = pList.get(5);
            int chromNum = pList.size()-6;
            this.chromLength = new int[chromNum];
            this.centStart = new int[chromNum];
            this.centEnd = new int[chromNum];
            for (int i = 0; i < chromNum; i++) {
                String[] tem = pList.get(i+6).split("\t");
                chromLength[i] = Integer.valueOf(tem[1]);
                centStart[i] = Integer.valueOf(tem[2]);
                centEnd[i] = Integer.valueOf(tem[3]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("Please select a parameter file with correct format");
            System.exit(1);
        }
    }
    
    public static void main (String[] args) {
        new DenovoEvaluation (args[0]);
        //new DenovoEvaluation (); //(test)
        
/* format of parameter file (args[0])       
        #working directory:
        /workdir/feilu/cml247Evaluation/
        #de novo assembly (Fasta or Fasta.gz format):
        /workdir/feilu/test.cml247.fasta.gz
        #genetic anchor file:
        /workdir/feilu/test.v2.panAnchor.rigid.togm.txt
        #minimum scaffold length (bp, default=20000, only scaffolds >= minimum scaffold length will be analyzed):
        20000
        #minimum anchor number on scaffold (default=20, only scaffolds with >= minimum anchor number will output)
        5
        #Reference genome information
        #Chromosome	Length(V2)	CentromereS	CentromereE
        1	301354134	134400000	135000000
        2	237068872	92900000	94700000
        3	232140173	99700000	100700000
        4	241473503	105300000	106100000
        5	217872851	102300000	109200000
        6	169174352	49600000	50000000
        7	176764761	54600000	62500000
        8	175793758	49000000	51400000
        9	156750705	72200000	72700000
        10	150189434	50100000	52400000
*/        
    }
}
