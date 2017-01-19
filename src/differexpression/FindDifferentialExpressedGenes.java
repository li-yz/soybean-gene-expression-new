package differexpression;

import utils.MyPrint;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by Liyanzhen on 2016/11/14.
 */
public class FindDifferentialExpressedGenes {
    public static void main(String [] args){
        List<String> GSEFileList = new ArrayList<>();
        GSEFileList.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE7108_Su_RMA_matrix.txt");
        GSEFileList.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE8432_Su_RMA_matrix.txt");//注意：：：： GSE8432这个系列实验按时间序列，现在选择的对照组与实验组可能不合适
        GSEFileList.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE29740_Su_RMA_matrix.txt");
        GSEFileList.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE29741_Su_RMA_matrix.txt");
        GSEFileList.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE33410_Su_RMA_matrix.txt");
        GSEFileList.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE41724_Su_RMA_matrix.txt");

        Set<String> differGenesOfAllGSEs = new HashSet<String>();

        FindDifferentialExpressedGenes find = new FindDifferentialExpressedGenes();
        for(String originSeriesDataPath :GSEFileList) {
            find.readMatrixDataAndFindDifferGenes(originSeriesDataPath, differGenesOfAllGSEs);
        }
        MyPrint.print("全部GSEs系列的差异表达基因",differGenesOfAllGSEs.size()+"");
        //保存找到的所有差异表达基因集合，把基因ID保存到txt文件中
        String outputPath = "D:\\paperdata\\soybean\\differExpressionGenes\\";
        outPutAllDifferGenesOfAllGSEs(outputPath,differGenesOfAllGSEs);

    }

    /**
     * 读取原始表达数据，保存到合适的数据结构里
     */
    public void readMatrixDataAndFindDifferGenes(String seriesMatrixDatapath ,Set<String> differGenesOfAllGSEs){
        //从文件路径中截取 实验系列名，类似GSE7108这样的名称
        String GSEName = "";
        Pattern pattern = Pattern.compile("(GSE){1}[0-9]*");
        Matcher matcher = pattern.matcher(seriesMatrixDatapath);
        while(matcher.find()){
            GSEName = matcher.group();
        }

        String[] genes = new String[61170];


        List<String> differGenesOfOneGSE = new ArrayList<String>();
        List<Integer> differIndexs = new ArrayList<Integer>();

        try {
            File file = new File(seriesMatrixDatapath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String [] all_rep_names = br.readLine().split("\t");
            int reps = all_rep_names.length;
            double [][] allGenesExpData = new double[61170][reps];

            //读取基因表达数据，genes保存所有基因探针ID、allGenesExpData保存读取到的基因表达数据
            SplitLine(genes, allGenesExpData, br);

            //找本次实验中差异表达基因集合
            findDifferGenes(GSEName,genes,allGenesExpData,differGenesOfOneGSE,differIndexs,differGenesOfAllGSEs);

            //文本文件 保存差异表达基因集合
            String outputPath = "D:\\paperdata\\soybean\\differExpressionGenes\\";
            outPutDifferGenesOfOneGSE(outputPath,differGenesOfOneGSE, GSEName);

        }catch (Exception e){
            System.out.println("读取基因表达数据文件异常");
            e.printStackTrace();
        }
    }

    /**
     * 读取原始的归一化之后的数据
     *
     * @param geneIds
     * @param br BufferdReader对象
     * @throws IOException
     */
    private void SplitLine(String[] geneIds, double[][] allGeneExpData, BufferedReader br) throws IOException {
        String line;
        int i = 0;
        while((line = br.readLine()) != null ){
            String [] oneLine = line.split("\t");
            geneIds[i] = oneLine[0];
            //
            for(int j=0; j < oneLine.length-1 ;j++){
                allGeneExpData[i][j] =Double.parseDouble(oneLine[j+1]);
            }

            i++;
        }
    }

    /**
     *计算每个基因对照组与实验组的fold-change值（比值版），f > 2的被选为差异表达基因
     *
     * @param genes 全部基因 ID
     * @param differGenes 用来保存找到的差异表达基因的ID
     * @param differIndexs 保存差异表达基因的下标
     */
    private void findDifferGenes(String GSEName,String[] genes,double[][] allGeneExpData,List<String>differGenes,List<Integer> differIndexs ,Set<String> differGenesOfAllGSEs){
        int gene_num = genes.length;
        double [] foldChange = new double[gene_num];



        for(int i=0;i < gene_num;i++){
            double meanValueOfControl;
            double meanValueOfTreatment;
            calculateFoldChangeForDifferGSE(GSEName , allGeneExpData, foldChange, i);
            if(isDifferExpressed(foldChange[i])){
                differGenes.add(genes[i]);
                differIndexs.add(i);
                differGenesOfAllGSEs.add(genes[i]);
            }
        }
        MyPrint.print(GSEName+"差异表达基因个数：","" + differGenes.size());
    }

    /**
     *判断一个基因是否 出现了差异表达。 可以灵活调整阈值，当前取2倍的 fold-change
     * @param v 一个基因i的fold change 值
     * @return 若v满足 差异表达的 阈值，则返回true, 否则返回 false
     */
    private boolean isDifferExpressed(double v) {
        return v >= 2 || v <= 0.5;
    }

    /**
     *
     * @param GSEName GSE系列名字，不同的系列的
     * @param allGeneExpData 归一化之后的表达数据，不同系列 的实验组/对照组 数据 分布规律不同，因此每个系列都要特定计算
     * @param foldChange 保存一个实验系列中 所有基因的fold-change值
     * @param i
     */
    private void calculateFoldChangeForDifferGSE(String GSEName, double[][] allGeneExpData, double[] foldChange, int i) {
        double meanValueOfControl =0.0d;
        double meanValueOfTreatment =0.0d;
        double sumControl = 0;
        double sumTreatment = 0;

        //由于不同的系列 GSE下，实验组、对照组数量不一、顺序不一，特定计算
        if("GSE7108".equals(GSEName)) {
            sumControl = allGeneExpData[i][0] + allGeneExpData[i][1] + allGeneExpData[i][4];
            meanValueOfControl = sumControl / 3;

            sumTreatment = allGeneExpData[i][2] + allGeneExpData[i][3] + allGeneExpData[i][5];
            meanValueOfTreatment = sumTreatment / 3;

            //由于RMA归一化之后的数据取了log,因此计算fold-change时，做 anti-log,
            foldChange[i] = Math.pow(2,meanValueOfControl-meanValueOfTreatment);
        }else if("GSE8432".equals(GSEName)){
            //GSE8432系列数据涉及到 2 个实验环境：Hw94-1 和 Tw72-1 ，1个 mock infection对照组
            //我使用此数据的方式： Hw94-1环境  VS  mock infection组  、  Tw72-1环境  VS  mock infection组，分别计算fold change值，找差异表达基因
            //注意：Hw94-1环境  与 Tw72-1环境各有一个芯片 取样失败

            double meanOfHw94_1 = 0.0d;
            double meanOfTw72_1 = 0.0d;
            double meanOfMock = 0.0d;

            meanOfMock = (allGeneExpData[i][0] + allGeneExpData[i][9] +allGeneExpData[i][18])/3;

            meanOfHw94_1 = (allGeneExpData[i][1] + allGeneExpData[i][3] + allGeneExpData[i][7] + allGeneExpData[i][10] + allGeneExpData[i][12]
                    + allGeneExpData[i][14] + allGeneExpData[i][16] + allGeneExpData[i][19] + allGeneExpData[i][21]+allGeneExpData[i][23] + allGeneExpData[i][25]) /11;
            meanOfTw72_1 = (allGeneExpData[i][2] + allGeneExpData[i][4] + allGeneExpData[i][8] + allGeneExpData[i][11] + allGeneExpData[i][13] + + allGeneExpData[i][15]
                    + allGeneExpData[i][17] + allGeneExpData[i][20] + + allGeneExpData[i][22] + allGeneExpData[i][24] + allGeneExpData[i][26]) /11;

            //由于RMA归一化之后的数据取了log,因此计算fold-change时，做 anti-log,
            double foldChangeHw94_1 = Math.pow(2,meanOfHw94_1 - meanOfMock);
            double foldChangeTw72_1 = Math.pow(2,meanOfTw72_1 - meanOfMock);

            //这两个 fold-change 值有一个在阈值范围内 就把该基因选作差异表达基因
            if(isDifferExpressed(foldChangeHw94_1)){
                foldChange[i] = foldChangeHw94_1;
            }else{
                foldChange[i] = foldChangeTw72_1;
            }

        }else if("GSE29740".equals(GSEName)){
            //GSE29740 实验设计：1种基因型的大豆：PI462312型；  3个实验条件：剧毒锈病环境、无毒性的锈病环境和mock infection ；  每一个环境下都有3个replications；  在感染后6个时间点取样
            //因此共有 1 * 3 * 3 * 6=54个基因芯片数据。
            //我使用此数据的方式：Tw90-2(vir)组  VS  mock infection  ，  Hw94-1(avir)组 VS mock infection 组



            double meanOfTw20 = 0.0d;
            double meanOfHw94 = 0.0d;
            double meanOfMockInfection=0.0d;

            double sumMock =0.0d;
            for(int j=0;j <= 17 ;j++){
                sumMock+=allGeneExpData[i][j];
            }
            meanOfMockInfection = (sumMock)/18;

            double sumTw20 = 0.0d;
            for(int j=18;j <=35;j++){
                sumTw20+=allGeneExpData[i][j];
            }
            meanOfTw20 = sumTw20/18;

            double sumHw94 =0.0d;
            for(int j=36;j <= 53;j++){
                sumHw94+=allGeneExpData[i][j];
            }
            meanOfHw94 = sumHw94/18;

            //由于RMA归一化之后的数据取了log,因此计算fold-change时，做 anti-log
            double foldChangeTw20 = Math.pow(2,meanOfTw20 -meanOfMockInfection);
            double foldChangeHw94 = Math.pow(2,meanOfHw94-meanOfMockInfection);

            //这两个 fold-change 值有一个在阈值范围内 就把该基因选作差异表达基因
            if(isDifferExpressed(foldChangeTw20)){
                foldChange[i] = foldChangeTw20;
            }else{
                foldChange[i] = foldChangeHw94;
            }

        }else if("GSE29741".equals(GSEName)){
            //GSE29741 的实验设计：2种基因型的大豆：PI459025B（有SBR Rpp4 resistant gene，即抗锈病） 和 Cultivar Williams 型（无任何SBR resistan gene，即不抗锈病）
            //2种 treatments：soybean rust isolate Hawaii 94-1 即含有锈病病菌的环境  和  mock infection （对照组）
            //每一种实验下都有3个replications
            //按时间序列采样，6个时间点，感染后12、24、72、144、216、288小时
            //因此总的 芯片数 = 2*2*3*6=72个
            //我的使用方式：不考虑时间因素， PI459025B 在 Hw94-1 和mock Infection计算一个 fold change, Cultivar Williams 型 也在这两个环境计算 一个fold change值

            double meanWilliamsMock = 0.0d;
            double meanWilliamsHw94_1 = 0.0d;
            double meanPI459025Mock = 0.0d;
            double meanPI459025Hw94_1 = 0.0d;

            double sumWilliamsMock = 0.0d;
            double sumWilliamsHw94_1 = 0.0d;
            double sumPI459025Mock = 0.0d;
            double sumPI459025Hw94_1 = 0.0d;

            for(int j=0;j <=17;j++){
                sumWilliamsMock+=allGeneExpData[i][j];
            }
            meanWilliamsMock = sumWilliamsMock/18;
            for(int j=18;j <=35;j++){
                sumWilliamsHw94_1+=allGeneExpData[i][j];
            }
            meanWilliamsHw94_1 =sumWilliamsHw94_1/18;

            for(int j=36;j <=53;j++){
                sumPI459025Mock+=allGeneExpData[i][j];
            }
            meanPI459025Mock = sumPI459025Mock/18;

            for(int j=54;j <=71;j++){
                sumPI459025Hw94_1 +=allGeneExpData[i][j];
            }
            meanPI459025Hw94_1 = sumPI459025Hw94_1/18;

            double foldChangePI459025Hw94_1 = Math.pow(2,meanPI459025Hw94_1-meanPI459025Mock);
            double foldChangeWilliams = Math.pow(2,meanWilliamsHw94_1-meanWilliamsMock);

            if(isDifferExpressed(foldChangePI459025Hw94_1)){
                foldChange[i] = foldChangePI459025Hw94_1;
            }else{
                foldChange[i] =foldChangeWilliams;
            }

        }else if("GSE33410".equals(GSEName)){
            //实验设计：2种基因型的大豆：PI23970（包含SBR Rpp2抗锈病基因） 和 Embraph型（Susceptible Brazilian cultivar 无抗锈病基因）
            //2treatments：soybean rust challenge（即锈病病菌环境） 和 mock infection （无锈病病菌环境 即对照组）
            //每个实验环境有3个replications
            //10个时间点 ，在感染后的6、12、18、24、36、48、72、96、120、168 小时后取样
            //因此，总的芯片个数 = 2*2*3*10 = 120个
            //我利用此数据的方式：①基因型 PI23970 锈病病菌环境  VS  mock infection   计算一个fold change
            //                    ②基因型 Embraph  锈病病菌环境  VS  mock infection   计算一个fold change

            double sumPI23970SBR =0.0d;
            double sumPI23970Mock =0.0d;
            double sumEmbraphSBR =0.0d;
            double sumEmbraphMock =0.0d;

            double meanPI23970SBR =0.0d;
            double meanPI23970Mock =0.0d;
            double meanEmbraphSBR =0.0d;
            double meanEmbraphMock =0.0d;

            for(int j=0;j <=56;j=j+3){
                sumPI23970SBR +=allGeneExpData[i][j++];
                sumPI23970SBR +=allGeneExpData[i][j++];
                sumPI23970SBR +=allGeneExpData[i][j++];
            }
            meanPI23970SBR = sumPI23970SBR/30;
            for(int j=3;j <=59;j=j+3){
                sumPI23970Mock += allGeneExpData[i][j++];
                sumPI23970Mock += allGeneExpData[i][j++];
                sumPI23970Mock += allGeneExpData[i][j++];
            }
            meanPI23970Mock = sumPI23970Mock/30;

            double foldChangePI23970 = Math.pow(2,meanPI23970SBR-meanPI23970Mock);

            for(int j=60;j <=116;j=j+3){
                sumEmbraphSBR += allGeneExpData[i][j++];
                sumEmbraphSBR += allGeneExpData[i][j++];
                sumEmbraphSBR += allGeneExpData[i][j++];
            }
            meanEmbraphSBR = sumEmbraphSBR/30;
            for(int j=63;j <=119;j=j+3){
                sumEmbraphMock += allGeneExpData[i][j++];
                sumEmbraphMock += allGeneExpData[i][j++];
                sumEmbraphMock += allGeneExpData[i][j++];
            }
            meanEmbraphMock = sumEmbraphMock/30;
            double foldChangeEmbraph = Math.pow(2,meanEmbraphSBR-meanEmbraphMock);
            if(isDifferExpressed(foldChangePI23970)){
                foldChange[i] = foldChangePI23970;
            }else{
                foldChange[i] = foldChangeEmbraph;
            }


        }else if("GSE41724".equals(GSEName)){
            //GSE41724 实验设计：次试验比较简单，只涉及到1种基因型的大豆，2个treatments
            double meanTreatment = (allGeneExpData[i][0]+allGeneExpData[i][1])/2;
            double meanControl = (allGeneExpData[i][2] + allGeneExpData[i][3] +allGeneExpData[i][4])/3;

            foldChange[i] = Math.pow(2,meanTreatment-meanControl);
        }

    }

    private void outPutDifferGenesOfOneGSE(String path, List<String> differGenesOfOneGSE, String GSEName){
        try {
            File file = new File(path+GSEName+"_differ_exp_genes.txt");
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));

            for(String geneId :differGenesOfOneGSE){
                bw.write(geneId+"\n");
            }
            bw.close();

        }catch (IOException e){
            e.printStackTrace();
        }
    }

    private static void outPutAllDifferGenesOfAllGSEs(String path, Set<String> differGenesOfAllGSEs){
        try {
            File file = new File(path+"all_differ_exp_genes.txt");
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));

            for(String geneId :differGenesOfAllGSEs){
                bw.write(geneId+"\n");
            }
            bw.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }

}
