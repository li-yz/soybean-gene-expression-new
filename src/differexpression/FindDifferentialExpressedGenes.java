package differexpression;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by Liyanzhen on 2016/11/14.
 */
public class FindDifferentialExpressedGenes {
    public static void main(String [] args){
        String originSeriesDataPath = "D:\\paperdata\\soybean\\originGSESeriesData\\GSE7108_series_matrix.txt";

        FindDifferentialExpressedGenes find = new FindDifferentialExpressedGenes();
        find.readSeriesMatrixData(originSeriesDataPath ,3 ,3);
    }

    /**
     * 读取原始表达数据，保存到合适的数据结构里
     */
    public void readSeriesMatrixData(String seriesMatrixDatapath,int controllRepNum, int treatmentRepNum){
        String[] genes = new String[61170];
        double[][] controlRep1icates = new double[61170][controllRepNum];
        double[][] treatmentReplicates = new double[61170][treatmentRepNum];
        List<String> differGenes = new ArrayList<String>();
        List<Integer> differIndexs = new ArrayList<Integer>();

        try {
            File file = new File(seriesMatrixDatapath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            while((line = br.readLine()) != null && !"!series_matrix_table_begin".equals(line)){

            }
            String [] table_title = br.readLine().split("\t");
            int i = 0;
            //读取基因表达数据，保存到genes、controlRep1icates、treatmentReplicates
            SplitLine(genes, controlRep1icates, treatmentReplicates, br, i);

            //找本次实验中差异表达基因集合
            findDifferGenes(genes,controlRep1icates,treatmentReplicates,differGenes,differIndexs);

            //文本文件 保存差异表达基因集合
            String path = "D:\\paperdata\\soybean\\differExpressionGenes\\GSE7108_differ_exp_genes.txt";
            outPutDifferGenes(path, differIndexs, genes, controlRep1icates, treatmentReplicates);

        }catch (Exception e){
            System.out.println("读取基因表达数据文件异常");
            e.printStackTrace();
        }
    }

    /**
     * 读取原始数据，按数据说明拆分数据（由于每个GSE序列的数据，对照组、实验组组数不一样，要针对每个文件重写这个拆分读取函数）
     * 适用于GSE7108的拆分方法
     * @param genes
     * @param controlRep1icates 保存对照组基因表达数据
     * @param treatmentReplicates 保存实验组基因表达数据
     * @param br BufferdReader对象
     * @param i 基因
     * @throws IOException
     */
    private void SplitLine(String[] genes, double[][] controlRep1icates, double[][] treatmentReplicates, BufferedReader br, int i) throws IOException {
        String line;
        while((line = br.readLine()) != null && !"!series_matrix_table_end".equals(line)){
            String [] geneExpression = line.split("\t");
            genes[i] = geneExpression[0];
            //对照组
            controlRep1icates[i][0] = Double.parseDouble(geneExpression[1]);
            controlRep1icates[i][1] = Double.parseDouble(geneExpression[2]);
            controlRep1icates[i][2] = Double.parseDouble(geneExpression[5]);

            //实验组
            treatmentReplicates[i][0] = Double.parseDouble(geneExpression[3]);
            treatmentReplicates[i][1] = Double.parseDouble(geneExpression[4]);
            treatmentReplicates[i][2] = Double.parseDouble(geneExpression[6]);
            i++;
        }
    }

    /**
     *计算每个基因对照组与实验组的fold-change值（比值版），f > 2的被选为差异表达基因
     * @param genes 全部基因 ID
     * @param controlRep1icates 对照组基因表达值数组
     * @param treatmentReplicates 实验组基因表达值数组
     * @param differGenes 用来保存找到的差异表达基因集合结果
     */
    private void findDifferGenes(String[] genes,double[][] controlRep1icates, double[][] treatmentReplicates,List<String>differGenes,List<Integer> differIndexs){
        int gene_num = genes.length;
        int controllRepNum=controlRep1icates[0].length;
        int treatmentRepNum=treatmentReplicates[0].length;
        double []meanValueOfControl = new double[gene_num];
        double []meanValueOfTreatment = new double[gene_num];
        double [] foldChange = new double[gene_num];

        for(int i=0;i < gene_num;i++){
            double sumControl = 0;
            for(int j=0;j < controllRepNum;j++){
                sumControl+=controlRep1icates[i][j];
            }
            meanValueOfControl[i] = sumControl/controllRepNum;
            double sumTreatment = 0;
            for(int j=0;j < treatmentRepNum;j++){
                sumTreatment+=treatmentReplicates[i][j];
            }
            meanValueOfTreatment[i] = sumTreatment/treatmentRepNum;
            foldChange[i] = meanValueOfControl[i]/meanValueOfTreatment[i];
            if(foldChange[i] >= 2 || foldChange[i] <= 0.5){
                differGenes.add(genes[i]);
                differIndexs.add(i);
            }
        }
    }

    private void outPutDifferGenes(String path,List<Integer> differIndexs,String[] genes,double [][]controlRep1icates,double [][]treatmentReplicates){
        try {
            File file = new File(path);
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));
            StringBuffer sb = new StringBuffer();
            sb.append("ID_REF");
            sb.append("\t");
            for(int t=0;t < treatmentReplicates[0].length;t++){
                sb.append("treatment");
                sb.append("\t");
            }
            for(int c=0;c < controlRep1icates[0].length;c++){
                sb.append("control");
                sb.append("\t");
            }
            bw.write(sb.toString() + "\n");
            for(int i=0;i < differIndexs.size();i++){
                int index = differIndexs.get(i);
                sb = new StringBuffer();
                sb.append(genes[index]);
                sb.append("\t");
                for(int t=0;t < treatmentReplicates[index].length;t++){
                    sb.append(treatmentReplicates[index][t]);
                    sb.append("\t");
                }
                for(int c=0;c < controlRep1icates[index].length;c++){
                    sb.append(controlRep1icates[index][c]);
                    sb.append("\t");
                }

                bw.write(sb.toString() + "\n");
            }


        }catch (IOException e){
            e.printStackTrace();
        }
    }

}
