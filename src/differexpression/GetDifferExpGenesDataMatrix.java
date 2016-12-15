package differexpression;

import utils.MyPrint;

import java.io.*;
import java.util.*;

/**
 * 前面FindDifferentialExpressedGenes 找到了差异表达基因 的基因ID集合。
 * 这一部分要做：利用前面的集合数据，和归一化之后的芯片数据，构造 差异表达基因数据 矩阵。行表示一个基因，列是每一个GSE实验条件的均值。
 *
 * Created by Liyanzhen on 2016/12/8.
 */
public class GetDifferExpGenesDataMatrix {
    public static void main(String[] args){
        GetDifferExpGenesDataMatrix obj = new GetDifferExpGenesDataMatrix();

        //先把找到的差异表达基因集合读取进来
        Set<String> differExpGenes = new HashSet<>();
        obj.readDifferGenes(differExpGenes);

        MyPrint.print("读取到差异表达基因ID集合","共："+differExpGenes.size());

        //开始构造 差异表达基因数据 矩阵。行表示一个基因，列是每一个GSE实验条件的均值(!!!暂时先取所有的replicates作为列！！！！)。
        //读取所有的 归一化之后的数据矩阵

        Map<String,List<Double>> differGenesExpData = new HashMap<>();//用来保存  “基因Id”---表达数据 数据表
        List<String> differExpDataTableHeader = new ArrayList<>();//用来保存数据表的 表头信息，即每一列 对应的是哪个实验条件，方便后续解释实验结果
        obj.startReadNormalizationData(differExpGenes,differGenesExpData,differExpDataTableHeader);
        MyPrint.print("构造差异表达基因数据 完成","行数："+differGenesExpData.size() +" 列数："+differExpDataTableHeader.size());


        //把构造出来的 差异表达基因数据 矩阵 保存到二进制文件中，用的时候可以直接读取到内存中（以二进制流形式，无需处理数据的格式）

    }

    public void startReadNormalizationData(Set<String> differExpGenes,Map<String,List<Double>> differGenesExpData, List<String> differExpDataTableHeader){
        List<String> dataPaths = new ArrayList<>();//所有GSE系列 归一化之后的数据文件 路径
        dataPaths.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE7108_Su_RMA_matrix.txt");
        dataPaths.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE8432_Su_RMA_matrix.txt");
        dataPaths.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE29740_Su_RMA_matrix.txt");
        dataPaths.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE29741_Su_RMA_matrix.txt");
        dataPaths.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE33410_Su_RMA_matrix.txt");
        dataPaths.add("D:\\paperdata\\soybean\\RMA normalization data\\GSE41724_Su_RMA_matrix.txt");
        for(String path:dataPaths){
            readAndConstructDifferExpData(differExpGenes,differGenesExpData,differExpDataTableHeader,path);
        }

    }
    private void readAndConstructDifferExpData(Set<String> differExpGenes,Map<String,List<Double>> differGenesExpData, List<String> differExpDataTableHeader,String path){
        try{
            File file = new File(path);
            BufferedReader br = new BufferedReader(new FileReader(file));

            //首先读取文件第一行，读取数据文件表头，保存
            for(String name: br.readLine().split("\t")){
                differExpDataTableHeader.add(name);
            }


            String line = "";
            while((line = br.readLine()) != null){
                String[] array = line.split("\t");
                if(differExpGenes.contains(array[0])) {//即如果当前基因是差异表达基因，才处理并保存该基因的表达数据
                    List<Double> geneExpData = null;
                    if (differGenesExpData.containsKey(array[0])) {
                        geneExpData = differGenesExpData.get(array[0]);
                    } else {
                        geneExpData = new ArrayList<>();
                    }
                    for (int i = 1; i < array.length; i++) {
                        geneExpData.add(Double.parseDouble(array[i]));
                    }

                    differGenesExpData.put(array[0], geneExpData);
                }
            }
            br.close();

        }catch (Exception e){
            MyPrint.print("读取数据文件"+path+" 异常");
            e.printStackTrace();
        }
    }



    public void readDifferGenes(Set<String> differExpGenes){
        String path ="D:\\paperdata\\soybean\\differExpressionGenes\\all_differ_exp_genes.txt";
        try{
            File file = new File(path);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            while((line = br.readLine()) != null){
                differExpGenes.add(line);
            }
            br.close();
        }catch (Exception e){
            MyPrint.print("读取之前找到的差异表达基因ID集合异常",path);
            e.printStackTrace();

        }
    }

}
