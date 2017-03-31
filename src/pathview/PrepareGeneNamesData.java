package pathview;

import differexpression.GetDifferExpGenesDataMatrix;
import utils.MyPrint;
import utils.MySerialization;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by Liyanzhen on 2017/3/29.
 */
public class PrepareGeneNamesData {
    public static void main(String[] args){
        GetDifferExpGenesDataMatrix getDifferExpGenesDataMatrix = new GetDifferExpGenesDataMatrix();
        //先把找到的差异表达基因集合读取进来
        Set<String> differExpGenes = new HashSet<>();
        getDifferExpGenesDataMatrix.readDifferGenes(differExpGenes);

        Map<String,String> geneIdMapEntrezId= getGeneIdMapEntrezId(differExpGenes);
        //序列化affyGeneId---entrez gene ID的映射关系
        MySerialization.serializeObject(geneIdMapEntrezId,"D:\\paperdata\\soybean\\基因ID-gene name\\geneIdMapEntrezId.obj");
        String path1 = "D:\\paperdata\\soybean\\基因ID-gene name\\all_DEGS_entrezId.txt";
        String header = "entrez";
        output(geneIdMapEntrezId,path1,header);

        Map<String,String> geneIdMapUNIGENE = getGeneIdMapUNIGENE(differExpGenes);
        String path2 = "D:\\paperdata\\soybean\\基因ID-gene name\\all_DEGS_UNIGENE.txt";
        header = "UNIGENE";
        output(geneIdMapUNIGENE,path2,header);

    }
    private static Map<String,String> getGeneIdMapEntrezId(Set<String> allDifferExpGenes){
        String path ="D:\\paperdata\\soybean\\基因ID-gene name\\geneIdToEntrezId.txt";
        Map<String,String> affyIdToEntrezId = new HashMap<>();
        try{
            File file = new File(path);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            String tableHeader = br.readLine();

            int count =0;
            while((line = br.readLine()) != null){
                String[]array = line.split("\t");
                if(array.length == 4){//说明DAVID的ID conversion工具成功的将原始的affy 类型的基因ID转换成了entrez_gene_id
                    String geneId = array[0];
                    String entrezId = array[1];
                    affyIdToEntrezId.put(geneId,entrezId);

                    if(allDifferExpGenes.contains(geneId)){
                        count++;
                    }
                }else {
                    //do nothing
                }
            }
            MyPrint.print("DAVID成功把原始文件中的Affimetrix 类型的基因ID转换成了Entrez Id,转换成功了："+count+"个");

            br.close();
        }catch (Exception e){
            MyPrint.print("读取geneId--转--entrezId文件异常",path);
            e.printStackTrace();

        }
        return affyIdToEntrezId;
    }

    private static Map<String,String> getGeneIdMapUNIGENE(Set<String> allDifferExpGenes){
        String path ="D:\\paperdata\\soybean\\基因ID-gene name\\affIdToUNIGENE.txt";
        Map<String,String> affyIdToEntrezId = new HashMap<>();
        try{
            File file = new File(path);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            String tableHeader = br.readLine();

            int count =0;
            while((line = br.readLine()) != null){
                String[]array = line.split("\t");
                if(array.length == 4){//说明DAVID的ID conversion工具成功的将原始的affy 类型的基因ID转换成了entrez_gene_id
                    String geneId = array[0];
                    String unigene = array[1];
                    affyIdToEntrezId.put(geneId,unigene);

                    if(allDifferExpGenes.contains(geneId)){
                        count++;
                    }
                }else {
                    //do nothing
                }
            }
            MyPrint.print("DAVID成功把原始文件中的Affimetrix 类型的基因ID转换成了UNIGENE,转换成功了："+count+"个");

            br.close();
        }catch (Exception e){
            MyPrint.print("读取geneId--转--entrezId文件异常",path);
            e.printStackTrace();

        }
        return affyIdToEntrezId;
    }

    private static void output(Map<String,String> geneIdMapName,String path,String header){
        try {
            FileWriter wr = new FileWriter(path);
            BufferedWriter br = new BufferedWriter(wr);

            br.write(header);
            br.newLine();
            for(Map.Entry entry :geneIdMapName.entrySet()){
                br.write((String)entry.getValue());
                br.newLine();
            }
            br.close();
            wr.close();
        }catch (Exception e){

        }
    }
}
