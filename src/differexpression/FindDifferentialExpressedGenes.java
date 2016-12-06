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

        Set<String> differGenesOfAllGSEs = new HashSet<String>();

        FindDifferentialExpressedGenes find = new FindDifferentialExpressedGenes();
        for(String originSeriesDataPath :GSEFileList) {
            find.readMatrixDataAndFindDifferGenes(originSeriesDataPath, differGenesOfAllGSEs);
        }
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
     * 读取原始数据，按数据说明拆分数据（由于每个GSE序列的数据，对照组、实验组组数不一样，要针对每个文件重写这个拆分读取函数）
     * 适用于GSE7108的拆分方法
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
    private void findDifferGenes(String GSEName,String[] genes,double[][] allGeneExpData,List<String>differGenes,List<Integer> differIndexs ,Set<String> differGenesOfAllGSMs){
        int gene_num = genes.length;
        double [] foldChange = new double[gene_num];



        for(int i=0;i < gene_num;i++){
            double meanValueOfControl;
            double meanValueOfTreatment;
            calculateFoldChangeForDifferGSE(GSEName , allGeneExpData, foldChange, i);
            if(foldChange[i] >= 2 || foldChange[i] <= 0.5){
                differGenes.add(genes[i]);
                differIndexs.add(i);
                differGenesOfAllGSMs.add(genes[i]);
            }
        }
        MyPrint.print("差异表达基因个数：","" + differGenes.size());
    }

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
        }else if("GSE8432".equals(GSEName)){
            sumControl = allGeneExpData[i][1] + allGeneExpData[i][3] + allGeneExpData[i][7] + allGeneExpData[i][10] + allGeneExpData[i][12]
                    + allGeneExpData[i][14] + allGeneExpData[i][16] + allGeneExpData[i][19] + allGeneExpData[i][21]+allGeneExpData[i][23] + allGeneExpData[i][25];
            meanValueOfControl = sumControl/11;
            sumTreatment = allGeneExpData[i][2] + allGeneExpData[i][4] + allGeneExpData[i][8] + allGeneExpData[i][11] + allGeneExpData[i][13] + + allGeneExpData[i][15]
                    + allGeneExpData[i][17] + allGeneExpData[i][20] + + allGeneExpData[i][22] + allGeneExpData[i][24] + allGeneExpData[i][26];
            meanValueOfTreatment = sumTreatment/11;
        }
        //由于RMA归一化之后的数据取了log,因此计算fold-change时，做 anti-log,
        foldChange[i] = Math.pow(2,meanValueOfControl-meanValueOfTreatment);
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
