package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by Liyanzhen on 2016/12/6.
 */
public class MathTest {
    public static void main(String [] args){
        double r = Math.pow(2,3);
        MyPrint.print("测试指数求法","" + r);



        Map<String,List<Integer>> differGenesExpData = new HashMap<>();

        List<Integer> list = new ArrayList<>();
        for(int i=0;i< 10;i++){
            list.add(i);
        }

        differGenesExpData.put("GmaAffx",list);

        System.out.println(differGenesExpData.get("GmaAffx"));
    }
}
