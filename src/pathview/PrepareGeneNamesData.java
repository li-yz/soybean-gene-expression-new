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
        //�Ȱ��ҵ��Ĳ�������򼯺϶�ȡ����
        Set<String> differExpGenes = new HashSet<>();
        getDifferExpGenesDataMatrix.readDifferGenes(differExpGenes);

        Map<String,String> geneIdMapEntrezId= getGeneIdMapEntrezId(differExpGenes);
        //���л�affyGeneId---entrez gene ID��ӳ���ϵ
        MySerialization.serializeObject(geneIdMapEntrezId,"D:\\paperdata\\soybean\\����ID-gene name\\geneIdMapEntrezId.obj");
        String path1 = "D:\\paperdata\\soybean\\����ID-gene name\\all_DEGS_entrezId.txt";
        String header = "entrez";
        output(geneIdMapEntrezId,path1,header);

        Map<String,String> geneIdMapUNIGENE = getGeneIdMapUNIGENE(differExpGenes);
        String path2 = "D:\\paperdata\\soybean\\����ID-gene name\\all_DEGS_UNIGENE.txt";
        header = "UNIGENE";
        output(geneIdMapUNIGENE,path2,header);

    }
    private static Map<String,String> getGeneIdMapEntrezId(Set<String> allDifferExpGenes){
        String path ="D:\\paperdata\\soybean\\����ID-gene name\\geneIdToEntrezId.txt";
        Map<String,String> affyIdToEntrezId = new HashMap<>();
        try{
            File file = new File(path);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            String tableHeader = br.readLine();

            int count =0;
            while((line = br.readLine()) != null){
                String[]array = line.split("\t");
                if(array.length == 4){//˵��DAVID��ID conversion���߳ɹ��Ľ�ԭʼ��affy ���͵Ļ���IDת������entrez_gene_id
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
            MyPrint.print("DAVID�ɹ���ԭʼ�ļ��е�Affimetrix ���͵Ļ���IDת������Entrez Id,ת���ɹ��ˣ�"+count+"��");

            br.close();
        }catch (Exception e){
            MyPrint.print("��ȡgeneId--ת--entrezId�ļ��쳣",path);
            e.printStackTrace();

        }
        return affyIdToEntrezId;
    }

    private static Map<String,String> getGeneIdMapUNIGENE(Set<String> allDifferExpGenes){
        String path ="D:\\paperdata\\soybean\\����ID-gene name\\affIdToUNIGENE.txt";
        Map<String,String> affyIdToEntrezId = new HashMap<>();
        try{
            File file = new File(path);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            String tableHeader = br.readLine();

            int count =0;
            while((line = br.readLine()) != null){
                String[]array = line.split("\t");
                if(array.length == 4){//˵��DAVID��ID conversion���߳ɹ��Ľ�ԭʼ��affy ���͵Ļ���IDת������entrez_gene_id
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
            MyPrint.print("DAVID�ɹ���ԭʼ�ļ��е�Affimetrix ���͵Ļ���IDת������UNIGENE,ת���ɹ��ˣ�"+count+"��");

            br.close();
        }catch (Exception e){
            MyPrint.print("��ȡgeneId--ת--entrezId�ļ��쳣",path);
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
